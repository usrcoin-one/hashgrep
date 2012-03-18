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


#define num_indexbits 20
#define num_tagbits 16
#define num_bfbits 128

typedef CuckooHashtable<uint64_t, uint32_t, num_tagbits, num_indexbits, num_bfbits> HashfilterType;
typedef HashfilterType::Bucket BucketType;

const uint32_t bucket_size = 4;
const uint32_t num_buckets = 1ULL << num_indexbits;
const uint32_t INDEXMASK = num_buckets - 1;
const uint32_t TAGMASK = (1ULL << num_tagbits) - 1;


HashfilterType hashfilter;

// texture<unsigned char, 1, cudaReadModeElementType> tex_buckets;
// texture<unsigned char, 1, cudaReadModeElementType> tex_strbuf;

unsigned char *dev_buckets;
unsigned char *dev_strbuf;
unsigned char *dev_corpus;
unsigned char *dev_pos_bitmap;

void exit_on_error(const char *name, cudaError_t err) {
    if (err) {
        if (err) printf("%s Error: %s\n", name, cudaGetErrorString(err));
        exit(-1);
    }
}

inline uint32_t  rol32(unsigned int word, int shift)
{
    return (word << shift) | (word >> (32 - shift));
}

inline __device__ uint32_t  dev_rol32(unsigned int word, int shift)
{
    return (word << shift) | (word >> (32 - shift));
}

__host__ void setup(char* pattern_file) 
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

    exit_on_error("cudaMalloc",
                cudaMalloc((void **) &dev_buckets, 
                           len));

    exit_on_error("cudaMalloc",
                cudaMalloc((void **) &dev_strbuf, 
                           strbuf_used));

    exit_on_error("cudaMemcpy",
                cudaMemcpy(dev_buckets, 
                           hashfilter.buckets_, 
                           len, 
                           cudaMemcpyHostToDevice));

    exit_on_error("cudaMemcpy",
                cudaMemcpy(dev_strbuf, 
                           strbuf, 
                           strbuf_used, 
                           cudaMemcpyHostToDevice));

 }


int pick_dim(dim3 &dimGrid,
            dim3 &dimBlock,
            int numthreads,
            int blocksize)
{
    unsigned int blocks_y = 1;
    unsigned int blocks_x = 1;
    unsigned int threads_1d = numthreads % blocksize;

    if (numthreads > (256 * blocksize)) {
        blocks_y = numthreads / (256 * blocksize);
        blocks_x = 256;
        threads_1d = blocksize;
    } else if (numthreads > blocksize) {
        blocks_x = numthreads / blocksize;
        threads_1d = blocksize;
    }

    unsigned int threads_used = blocks_y * blocks_x * threads_1d;
    numthreads -= threads_used;
    printf("pick_dim %d %d %d\n", blocks_y, blocks_x, threads_1d);
    dimGrid = dim3(blocks_x, blocks_y);
    dimBlock = dim3(threads_1d);
    return threads_used;
}


inline __device__ void dev_set_bit(int i, unsigned char *bv) {
    unsigned int *p = (unsigned int *)bv;
    unsigned int bitMask = 1 << (i & 31);
    atomicOr(&p[i >> 5], bitMask);
}

inline __device__ uint32_t dev_read_tag(unsigned char* dev_buckets,
                                        uint32_t i, 
                                        uint32_t j)
{
    size_t offset = num_tagbits * j;
    HashfilterType::Bucket* b;
    b = (HashfilterType::Bucket*) (dev_buckets + i * sizeof(HashfilterType::Bucket) );
    uint32_t v = *(uint32_t *) (b->tagbits_ + offset / 8);
    v = (v >> (offset & 0x7)) & TAGMASK;
    return v;
}

inline __device__ uint32_t dev_read_value(unsigned char* dev_buckets, 
                                          uint32_t i, 
                                          uint32_t j)
{
    HashfilterType::Bucket* b;
    b = (HashfilterType::Bucket*) (dev_buckets + i * sizeof(HashfilterType::Bucket) );
    return b->valbits_[j];
}


__global__ void GrepKernel(unsigned char *d_a,
                           unsigned char *dev_buckets,
                           unsigned char *dev_pos_bitmap,
                           unsigned int char_offset)
{
    __shared__ unsigned boxed[BLOCK_SIZE + HASH_LEN];

    int i = char_offset + (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x * blockDim.y + threadIdx.x;

    boxed[threadIdx.x] = sbox[d_a[i]];

    /* Ugly, but let some threads pull in the remainder */
    /* TIME:  0.01 seconds */
    int otid = threadIdx.x;
    if (otid < HASH_LEN) {
        int new_i = blockDim.x + i;
        int new_b = blockDim.x + otid;
        boxed[new_b] = sbox[d_a[new_i]];
    }

    /* TIME:  Almost none.  */
    __syncthreads();

    uint32_t  hval[2] = {0, 0};
    /* Step 2:  Compute the hash of the next HASH_LEN characters */
    for (int j = 0; j < HASH_LEN; j++) {
        hval[0] = dev_rol32(hval[0], 1);
        hval[1] = dev_rol32(hval[1], 3);
        uint32_t sbv = boxed[threadIdx.x+j];
        hval[0] ^= sbv;
        hval[1] ^= sbv;
    }


    uint32_t i1, i2, tag;
    bool found = false;
    
    tag =  hval[0] & TAGMASK;
    tag += (tag == 0); 
    i1 = hval[1] & INDEXMASK;
    i2 = (i1 ^ (tag * 0x5bd1e995)) & INDEXMASK;

    BucketType *b1 = (BucketType *) (dev_buckets + i1 * sizeof(BucketType));

    

    uint64_t tagbits1 = *(uint64_t*) b1->tagbits_;

    for (int j = 0; j < bucket_size; j ++) {
        
        uint32_t tag1 = tagbits1 & TAGMASK;
        tagbits1 >>= num_tagbits;
        //char tag1 = b1->tagbits_[j];
        if (tag1 == tag) {
            found = true;
            break;
        }

        // BucketType *b2 = (BucketType *) (dev_buckets + i2 * sizeof(BucketType));
        // uint32_t tag2 = *((uint32_t*) b2->tagbits_);
        // if (tag2 == tag) {
        //     found = true;
        //     break;
        // }


        //uint32_t val;
        // if (dev_read_tag(dev_buckets, i1, j) == tag) {
        //     //val = dev_read_value(i1, j);
        //     found = true;
        //     break;
        // }

        // if (dev_read_tag(dev_buckets, i2, j) == tag) {
        //     //val = dev_read_value(i2, j);
        //     found = true;
        //     break;
        // }        
    }

    if (found) {
        dev_set_bit(i, dev_pos_bitmap);
    }
}

__host__ void process_corpus(char* corpus_file) 
{
    char *pinnedBuf;

    printf("opening corpus file %s\n", corpus_file);
    size_t filesize;
    int f = getfile(corpus_file, &filesize);
    if (f == -1) {
        perror(corpus_file);
        exit(-1);
    }

    filesize = min((unsigned long long)filesize, (unsigned long long)FILE_MAX);

    exit_on_error("cudaMallocHost pinnedBuf", 
                  cudaMallocHost((void **)&pinnedBuf, filesize));

    exit_on_error("cudaMalloc dev_corpus",
                  cudaMalloc((void **)&dev_corpus, filesize + HASH_LEN));

    exit_on_error("cudaMalloc dev_pos_bitmap",
                  cudaMalloc((void **)&dev_pos_bitmap, filesize / 8 + 1));

    exit_on_error("cudaMemset dev_pos_bitmap = 0",
                  cudaMemset(dev_pos_bitmap, 0, filesize/8 + 1));


    int numthreads;
    dim3 dimGrid, dimBlock;

    cudaStream_t streams[NR_STREAMS];

    for (int i = 0; i < NR_STREAMS; i++) {
        exit_on_error("cudaStreamCreate",
                    cudaStreamCreate(&streams[i]));
    }

    int size = filesize / NR_STREAMS;

    int fd = open(corpus_file, O_RDONLY);

    for (int i = 0; i < NR_STREAMS; i++) {
        unsigned offset = i * size;
        if (i == NR_STREAMS - 1) {
            size = filesize - i * size;
        }
        numthreads = size;

        printf("Executing grep on %d\n", size);

        read(fd, pinnedBuf + offset, size);

        exit_on_error("cudaMemcpyAsync",
            cudaMemcpyAsync(dev_corpus + offset, pinnedBuf + offset, size, cudaMemcpyHostToDevice, streams[i]));

        unsigned int char_offset = 0;
        while (numthreads > 0) {
            unsigned int tu = pick_dim(dimGrid, dimBlock, numthreads, BLOCK_SIZE);
            printf("Executing GrepKernel (%d,%d,%d) @ %u\n", dimGrid.x, dimGrid.y, dimBlock.x, (offset + char_offset));
            GrepKernel<<<dimGrid, dimBlock, 0, streams[i]>>>(dev_corpus, 
                                                             dev_buckets,
                                                             dev_pos_bitmap, 
                                                             offset + char_offset);
            //checkReportCudaStatus("GrepKernel");
            numthreads -= tu;
            char_offset += tu;
        }
    }

    cudaThreadSynchronize();

    for (int i = 0; i < NR_STREAMS; i++) {
        cudaStreamDestroy(streams[i]);
    }
    close(f);
    close(fd);

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

    setup(patterns_file);

    timing_stamp("setup", false);

    process_corpus(corpus_file);

    timing_stamp("processCorpus", false);

    timing_stamp("cleanup done", true);

    timing_report();

    struct cudaDeviceProp cdp;
    cudaGetDeviceProperties(&cdp, 0);
    printf("\ndeviceOverlap = %d\n", cdp.deviceOverlap);
}

