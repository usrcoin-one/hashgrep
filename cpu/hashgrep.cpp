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

Filter::Filter() {}

Filter::~Filter() {}


Filter::updateHashes(u_int8_t nextChar)
{
    if (nextChar & 0x80)
    {
        nextChar = 0;
    }

    hash.update(nextChar);
}

Filter::operator<<(string& pattern) throw(FilterError)
{
    if (pattern.length() < DEF_PATT_LEN) {
        return *this;
        //throw PatternError("Search phrase too short");
    }

    hash.reset();

    for (unsigned i = 0; i < DEF_PATT_LEN; i++) {
        updateHashes(pattern[i]);
    }

    uint32_t kbuf[2];
    kbuf[0] = hash.hval1();
    kbuf[1] = hash.hval2();
    Value key( (char *)kbuf, sizeof(kbuf) );
    hashfilter.Add(key);
    assert(hashfilter.Contain(key) == Ok);
    return *this;
}

void Filter::processCorpus(int fd) throw(FilterError)
{
    const int size = 512 * 1024;
    char* buff = new char[size];

    hash.reset();

    while (true) {
        int r = read(fd, (void*)buff, size);
        if (r < 0) {
            perror("Error");
            throw FilterError("Error while reading from file");
        }

        if (r == 0) {
            break;
        }

        for (int i = 0; i < r; i++) {

            if (buff[i] == '\n') {
                hash.reset();
                continue;
            }

            updateHashes(buff[i]);

            if (hash.is_full()) {

                if (hashfilter.Get(hash.h, ) == Ok) {
                    
                }
            }
        }

    }

    delete [] buff;
}

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cerr << "Usage: "<< argv[0] <<" phrases < corpus > filtered_corpus" << endl;
        return -1;
    }


    Filter<hash_rot_sbox_pre_2<DEF_PATT_LEN>, DEF_PATT_LEN> filter;

    //
    // Generate hashfilter from phrases
    //
    ifstream inPhrases(argv[1]);
    string phrase;
    while (getline(inPhrases, phrase))
    {
        try
        {
            filter << phrase;
        }
        catch (FilterError& fe)
        {
            cerr << "------at phrase:" << endl;
            cerr << phrase << endl;
            cerr << "Exception " << fe.what() << endl;
        }
    }
    inPhrases.close();

    filter.printSummary(STDOUT_FILENO);

    try
    {
        TIME("processFile", filter.processCorpus(STDIN_FILENO));
        //TIME("filterPhrases", filter.filterPhrases(pfd, outPhrases, outPhrasesIndex));
    }
    catch (FilterError& fe)
    {
        cerr << "Exception " << fe.what() << endl;
    }

    //cerr << "fprate: " << contor << endl;

    return 0;
}
