/**
Copyright 2012 Carnegie Mellon University

Authors: Bin Fan, Iulian Moraru and David G. Andersen

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This is the implementation of a feed-forward Bloom filter.
The code has been tested on a Linux 2.6.35 machine, gcc version 4.4.5.
*/

#include <iostream>
#include <stdio.h>
#include <errno.h>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <xmmintrin.h>
#include <random>

#include "hashgrep.h"


#define DO_TIMING 1
#if DO_TIMING
#define TIME(label, statement) \
    do { \
    struct timeval tvs, tve; \
    gettimeofday(&tvs, NULL); \
    do { statement; } while(0);	\
    gettimeofday(&tve, NULL); \
    double tvsd = (double)tvs.tv_sec + (double)tvs.tv_usec/1000000; \
    double tved = (double)tve.tv_sec + (double)tve.tv_usec/1000000; \
    fprintf(stderr, "%s time: %.5f\n", label, tved-tvsd); \
    } while (0)
#else
#define TIME(label, statement) statement
#endif


#ifndef DEF_CACHE_SIZE
#define DEF_CACHE_SIZE 0x1000000
#endif

#ifndef HUGEPAGE_FILE_NAME
#define HUGEPAGE_FILE_NAME "./mnt/hugepagefile"
#endif

const uint32_t Filter::BloomCacheMask = DEF_CACHE_SIZE - 1;

using namespace std;

mt19937 engine;

Filter::Filter() {
    //bitvector = new uint8_t[DEF_CACHE_SIZE];
    //memset(bitvector, 0, DEF_CACHE_SIZE);
    //LARGE PAGE
    fid = open(HUGEPAGE_FILE_NAME, O_CREAT | O_RDWR, 0755);
    if (fid < 0)
    {
		cerr << HUGEPAGE_FILE_NAME << endl;
        perror("open file error: ");
        throw FilterError("Not able to initialize filter");
    }

    mapaddr = mmap(NULL, DEF_CACHE_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fid, 0);

    if (mapaddr == MAP_FAILED)
    {
        perror("map failed: ");
        throw FilterError("Not able to initialize filter");
    }

    bitvector = (uint8_t*)mapaddr;
}

Filter::~Filter() {
    //delete [] bitvector;    
    munmap(mapaddr, DEF_CACHE_SIZE);
    close(fid);
    unlink(HUGEPAGE_FILE_NAME);
}


inline void Filter::updateHashes(u_int8_t nextChar)
{
    if (nextChar & 0x80)
    {
        nextChar = 0;
    }

    hash.update(nextChar);
}

inline void Filter::setBit(uint32_t index)
{
    uint32_t byteIndex = index >> 3;
    uint8_t bitMask = 1 << (index & 0x00000007);
    bitvector[byteIndex] |= bitMask;
}

inline bool Filter::checkInFilter()
{
    uint32_t x1 = hash.hval1() & BloomCacheMask;
    uint32_t x2 = hash.hval2() & BloomCacheMask;
    uint32_t byteIndex;
    uint8_t bitMask;

    byteIndex= x1 >> 3;
    bitMask = 1 << (x1 & 0x00000007);
    if (!(bitvector[byteIndex] & bitMask))
    {
        return false;
    }

    byteIndex = x2 >> 3;
    bitMask = 1 << (x2 & 0x00000007);
    if (!(bitvector[byteIndex] & bitMask))
    {
        return false;
    }

    return true;
}

void Filter::buildFilter(char* phrasefilename) throw(FilterError)
{
    ifstream inPhrases(phrasefilename);
    if (!inPhrases.is_open()) {
        cerr << "Can not open phrase file " << phrasefilename << endl;
        exit(-1);
    }

    string pattern;
    while (getline(inPhrases, pattern))
    {
        try
        {

            if (pattern.length() < DEF_PATT_LEN) {
                throw FilterError("Search phrase too short");
            }

            hash.reset();

            for (unsigned i = 0; i < DEF_PATT_LEN; i++) {
                updateHashes(pattern[i]);
            }

            char* p;
            if (hashfilter.Get(hash.h, p) == Ok) {
        
            } else
                p = NULL;

            char* cstr = new char [pattern.size() + sizeof(char*) + 1];
            memcpy(cstr, &p, sizeof(char*));
            strcpy(cstr + sizeof(char*), pattern.c_str());

            if (hashfilter.Put(hash.h, cstr) != Ok) {
                throw FilterError("Error while buiding hash table");
            }

            uint32_t x1 = hash.hval1() & BloomCacheMask;
            uint32_t x2 = hash.hval2() & BloomCacheMask;
            setBit(x1);
            setBit(x2);

        }
        catch (FilterError& fe)
        {
            cerr << "------at phrase:" << endl;
            cerr << pattern << endl;
            cerr << "Exception " << fe.what() << endl;
            break;
        }
    }
    inPhrases.close();
}

void Filter::processCorpus(int fd) throw(FilterError)
{
    const int size = 512 * 1024;
    char* buff = new char[size];

    string line;
    size_t total = 0;
    size_t hit = 0;
    while (getline(cin, line))
    {

        //cout << "----" << line << endl;
        bool printLine = false;
        hash.reset();
        const char* buff = line.data();
        char* hashstart = (char* ) buff;
        for (unsigned i = 0; i < line.length(); i++) {
            updateHashes(buff[i]);
            
            if (hash.is_full()) {
                //cout << "\t with " << i << "  \"" << hashstart << "\"" << endl;
                total ++;
                if (!checkInFilter()) {
                    hashstart ++;
                    continue;
                }

                hit ++; 
                
                //hashstart ++; continue;

                char* p;
                if (hashfilter.Get(hash.h, p) == Ok) {

                    while (p != NULL) {
                        char* q = p + sizeof(char*);
                        //cout << "\t checking " << q << endl;
                        if (strncmp(hashstart, q, strlen(q)) == 0) {
                            printLine = true;
                            break;
                        }
                        p = *(char**) p;
                    }

                    if (printLine) {
                        cout << line << endl;
                        break;
                    }
                }
                hashstart  ++;
               
            }
        }

    }

    cout << total << " " << hit << endl;

    delete [] buff;
}

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cerr << "Usage: "<< argv[0] <<" phrases < corpus > filtered_corpus" << endl;
        return -1;
    }


    Filter filter;



    try
    {
        TIME("buildFilter", filter.buildFilter(argv[1]));
        filter.printStatus();
        TIME("processCorpus", filter.processCorpus(STDIN_FILENO));
    }
    catch (FilterError& fe)
    {
        cerr << "Exception " << fe.what() << endl;
    }

    return 0;
}
