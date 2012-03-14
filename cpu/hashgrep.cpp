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

This is the implementation of a hash based grep
The code has been tested on a Linux 2.6.38 machine, gcc version 4.6.1.
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
#include <vector>

#include "hashgrep.h"

#define HUGEPAGE

#ifndef DEF_CACHE_SIZE
#define DEF_CACHE_SIZE 0x1000000
#endif

#ifndef HUGEPAGE_FILE_NAME
#define HUGEPAGE_FILE_NAME "./mnt/hugepagefile"
#endif

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


const uint32_t Filter::BloomCacheMask = DEF_CACHE_SIZE - 1;

using namespace std;

Filter::Filter() {

#ifdef HUGEPAGE
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
#else
    bitvector = new uint8_t[DEF_CACHE_SIZE];
    memset(bitvector, 0, DEF_CACHE_SIZE);
#endif

    strbuf_size = 0;
}

Filter::~Filter() {
    
#ifdef HUGEPAGE 
    munmap(mapaddr, DEF_CACHE_SIZE);
    close(fid);
    unlink(HUGEPAGE_FILE_NAME);
#else
    delete [] bitvector;   
#endif
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

void Filter::processPatterns(char* phrasefilename) throw(FilterError)
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
            strbuf_size += pattern.size() + sizeof(char*) + 1;

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

    int crtLineSize = 4*1024;
    char* line = new char[crtLineSize + 1];
    line[0] = '\0';

    size_t cnt0 = 0, cnt1 = 0, cnt2 = 0, cnt3 = 0;

    vector<int> vec1;
    vector<char*> vec2;
    int linePos = 0;

    hash.reset();
    while(true) {
        int r = read(fd, (void*)buff, size);
        if (r < 0) {
            perror("Error");
            throw FilterError("Error while reading from file");
        }
        if (r == 0) {
            break;
        }

        int lineStart = 0;

        for (int i = 0; i < r; i++) {
            if (buff[i] == '\n') {
                buff[i] = '\0';
                char* start;
                if (line[0] != '\0') {
                    // line is in use:
                    strncat(line, buff + lineStart, i - lineStart + 1 );
                    start = line;
                }
                else {
                    start = buff + lineStart;
                }

                // perform strcmp for all recorded
                bool printLine = false;
                for (unsigned j = 0; j < vec1.size(); j ++ ) {
                    char *pos = start + vec1[j] - DEF_PATT_LEN + 1;
                    char *p = vec2[j];
                    while (p != NULL) {
                        char* q = p + sizeof(char*);
                        if (strncmp(pos, q, strlen(q)) == 0) {
                            printLine = true;
                            break;
                        }
                        p = *(char**) p;
                    }
                    if (printLine) {
                        cnt3 ++;
                        cout << start << endl;
                        break;
                    }
                    
                }
                vec1.clear();
                vec2.clear();
                line[0] = '\0';

                hash.reset();
                lineStart = i + 1;
                linePos = 0;
                continue;
            }
          
            updateHashes(buff[i]);
            //cout << buff[i] << " i = " << i << " linePos = " << linePos << " lineStar = " << lineStart << endl;

            if (hash.is_full()) {
                cnt0 ++;
                if (checkInFilter()) {
                    cnt1 ++;

                    char* p;
                    if (hashfilter.Get(hash.h, p) == Ok) {
                        cnt2 ++;
                        vec1.push_back(linePos);
                        vec2.push_back(p);                        
                    } 
                }
            }
            
            linePos ++;
        }

        if (buff[r - 1] != '\n') {
            if (r - lineStart > crtLineSize) {
                crtLineSize = (r - lineStart) * 2;
                delete [] line;
                line = new char[crtLineSize + 1];
            }
            memcpy(line, buff + lineStart, r - lineStart);
            line[r - lineStart] = '\0';
        }

    }

    cerr << "Total " << cnt0 << " queries" << endl;
    cerr << "Bloomfilter sees " << cnt1 << " hits, hitratio=" << 1.0 * cnt1/cnt0 << endl;
    cerr << "Hashfilter sees " << cnt2 << " hits, hitratio=" << 1.0 * cnt2/cnt1 << endl;
    cerr << "Strcmp sees " << cnt3 << " hits" << endl;
    delete [] buff;
}


void Filter::printStatus() {
    cerr << "Bloomfilter size  =" << (DEF_CACHE_SIZE >> 10) << " KB" << endl;
    cerr << "Hashfilter size   =" << (hashfilter.SizeInBytes() >> 10) << " KB" << endl;
    cerr << "All string buffer =" << (strbuf_size >> 10) << " KB" << endl;
}

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cerr << "Usage: "<< argv[0] <<" phrases < corpus > filtered_corpus" << endl;
        return -1;
    }


    Filter filter;

    try {
        TIME("processPatterns", filter.processPatterns(argv[1]));
        TIME("processCorpus", filter.processCorpus(STDIN_FILENO));
        filter.printStatus();
    }
    catch (FilterError& fe) {
        cerr << "Exception " << fe.what() << endl;
    }

    return 0;
}
