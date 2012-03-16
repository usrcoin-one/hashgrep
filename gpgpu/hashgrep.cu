/* -*- Mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*- */

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

This is the implementation of a hash based grep for GPGPU.
*/

#define _DARWIN_FEATURE_64_BIT_INODE 1
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
 #include <unistd.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>

#include "cuda.h"
#include "cuda_runtime.h"
#include "math_functions.h"

#include "sbox.h"
#include "cuckoohashtable.h"

using namespace std;


extern "C" {
int getfile(char *infile, size_t *filesize);
#include "timing.h"
}

#define BLOCK_SIZE 256
#define HASH_LEN 19
#define FILE_MAX 6710886400
#define NR_STREAMS 10

#ifndef DEF_STRBUF_SIZE
#define DEF_STRBUF_SIZE 0xa000000
#endif


char *pinnedBuf;
char *strbuf;
size_t strbuf_used = 1;

typedef CuckooHashtable<uint64_t, uint32_t, 16, 20, 128> HashfilterType;

HashfilterType hashfilter;

texture<unsigned char, 1, cudaReadModeElementType> tex_buckets;
texture<unsigned char, 1, cudaReadModeElementType> tex_strbuf;
unsigned char *dev_buckets;
unsigned char *dev_strbuf;

void exitOnError(const char *name, cudaError_t err) {
    if (err) {
        if (err) printf("%s Error: %s\n", name, cudaGetErrorString(err));
        exit(-1);
    }
}

inline uint32_t  rol32(unsigned int word, int shift)
{
    return (word << shift) | (word >> (32 - shift));
}

__host__ void process_patterns(char* pattern_file) 
{
    ifstream inPhrases(pattern_file);
    if (!inPhrases.is_open()) {
        cerr << "Can not open phrase file " << pattern_file << endl;
        exit(-1);
    }

    string pattern;
    while (getline(inPhrases, pattern)) {
        if (pattern.length() < HASH_LEN) {
            perror("Search phrase too short");
            exit(-1);
        }

        uint32_t  hval[2] = {0, 0};

        for (int j = 0; j < HASH_LEN; j++) {
            hval[0] = rol32(hval[0], 1);
            hval[1] = rol32(hval[1], 3);
            uint32_t  sbv = sbox[pattern[j]];
            hval[0] ^= sbv;
            hval[1] ^= sbv;
        }
            
        uint32_t p;
        if (hashfilter.Get(*((uint64_t*) hval), p) == Ok) {
            assert(p > 0);
        } else
            p = 0;

        if (strbuf_used >= DEF_STRBUF_SIZE)  {
             perror("Not enough strbuf, please make DEF_STRBUF_SIZE larger ");
        }

        char* cstr = strbuf + strbuf_used;
        memcpy(cstr, &p, sizeof(uint32_t));
        strcpy(cstr + sizeof(uint32_t), pattern.c_str());

        if (hashfilter.Put(*((uint64_t*) hval), strbuf_used) != Ok) {
            perror("Error while buiding hash table");
            exit(-1);
        }
        strbuf_used += pattern.size() + sizeof(uint32_t) + 1;
            
    }
    inPhrases.close();
   
    hashfilter.BuildBF();

    size_t len = hashfilter.SizeInBytes();

    exitOnError("cudaMalloc",
                cudaMalloc((void **) &dev_buckets, 
                           len));

    exitOnError("cudaMalloc",
                cudaMalloc((void **) &dev_strbuf, 
                           strbuf_used));

    exitOnError("cudaMemcpy",
                cudaMemcpy(dev_buckets, 
                           hashfilter.buckets_, 
                           len, 
                           cudaMemcpyHostToDevice));

    exitOnError("cudaMemcpy",
                cudaMemcpy(dev_strbuf, 
                           strbuf, 
                           strbuf_used, 
                           cudaMemcpyHostToDevice));

 }


__host__ void process_corpus(char* corpus_file) 
{

    printf("GPUGrep opening %s\n", corpus_file);
    size_t filesize;
    int f = getfile(corpus_file, &filesize);
    if (f == -1) {
        perror(corpus_file);
        exit(-1);
    }
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        cerr << "usage: " << argv[0] << " patterns corpus" << endl;
        return -1;
    }

    char *patterns_file = argv[1];
    char *corpus_file = argv[2];

    //cout << sizeof(HashfilterType::Bucket) << endl;

    strbuf = new char[DEF_STRBUF_SIZE];

    timing_stamp("start", false);

    process_patterns(patterns_file);

    timing_stamp("setup", false);

    process_corpus(corpus_file);

    timing_stamp("processCorpus", false);

    timing_stamp("cleanup done", true);

    timing_report();

    struct cudaDeviceProp cdp;
    cudaGetDeviceProperties(&cdp, 0);
    printf("\ndeviceOverlap = %d\n", cdp.deviceOverlap);
}

