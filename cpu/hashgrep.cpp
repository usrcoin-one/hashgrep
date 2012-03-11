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


using namespace std;

mt19937 engine;

Filter::Filter() {}

Filter::~Filter() {}


void Filter::updateHashes(u_int8_t nextChar)
{
    if (nextChar & 0x80)
    {
        nextChar = 0;
    }

    hash.update(nextChar);
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
                //cout << "\t with " << i << "  " << hashstart << endl;

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
                hashstart ++;
            }
        }

    }
    //while (true) {
        // int r = read(fd, (void*)buff, size);
        // if (r < 0) {
        //     perror("Error");
        //     throw FilterError("Error while reading from file");
        // }

        // if (r == 0) {
        //     break;
        // }

        // for (int i = 0; i < r; i++) {

        //     if (buff[i] == '\n') {
        //         hash.reset();
        //         continue;
        //     }

        //     updateHashes(buff[i]);

        //     if (hash.is_full()) {
        //         char* p;
        //         if (hashfilter.Get(hash.h, p) == Ok) {
        //             //cout << "have seen this prefix:" << endl;
        //             patterns.push_back(p);
        //             linestarts.push_back(
    
        //     }
        // }

    //}

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
        TIME("processFile", filter.processCorpus(STDIN_FILENO));
    }
    catch (FilterError& fe)
    {
        cerr << "Exception " << fe.what() << endl;
    }

    return 0;
}
