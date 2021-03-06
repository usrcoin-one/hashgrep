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
//#define NR_STREAMS 10
#define NR_STREAMS 1

#ifndef DEF_STRBUF_SIZE
#define DEF_STRBUF_SIZE 0xa000000
#endif

//#define ENABLE_COUNTING
#define ENABLE_FILTERING
#define MEM_ALIGNMENT 8


char *pinnedBuf;
char *strbuf;
uint32_t strbuf_used = MEM_ALIGNMENT;


#define num_indexbits 20
#define num_tagbits 16
#define num_bfbits 128

typedef CuckooHashtable<uint64_t, uint32_t, num_tagbits, num_indexbits, num_bfbits> HashfilterType;
typedef HashfilterType::Bucket BucketType;

const uint32_t bucket_size = 4;
const uint32_t num_buckets = 1ULL << num_indexbits;
const uint32_t INDEXMASK = num_buckets - 1;
const uint32_t TAGMASK = (1ULL << num_tagbits) - 1;
const uint32_t BFINDEXMASK  = num_bfbits - 1;

HashfilterType hashfilter;

unsigned char *dev_buckets;
unsigned char *dev_strbuf;
unsigned char *dev_corpus;
unsigned char *dev_pos_bitmap;
unsigned int *dev_cnt1;
unsigned int *dev_cnt2;

void exit_on_error(const char *name, cudaError_t err) {
    if (err) {
        if (err) printf("%s Error: %s\n", name, cudaGetErrorString(err));
        exit(-1);
    }
}


void checkReportCudaStatus(const char *name) {
    cudaError_t err = cudaGetLastError();
    printf("CudaStatus %s: ", name);
    if (err) printf("Error: %s\n", cudaGetErrorString(err));
    else printf("Success\n");
}


inline uint32_t  rol32(unsigned int word, int shift)
{
    return (word << shift) | (word >> (32 - shift));
}

inline __device__ uint32_t  dev_rol32(unsigned int word, int shift)
{
    return (word << shift) | (word >> (32 - shift));
}

void setup(char* pattern_file) 
{
    ifstream inPhrases(pattern_file);
    if (!inPhrases.is_open()) {
        cerr << "Can not open phrase file " << pattern_file << endl;
        exit(-1);
    }

    char* logbuf = new char[DEF_STRBUF_SIZE];
    uint32_t logbuf_used = MEM_ALIGNMENT;
    if (!logbuf) {
        perror("cannot allocate host memory for strbuf\n");
        exit(-1);
    }
    memset(logbuf, 0, DEF_STRBUF_SIZE);

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
            uint32_t  sbv = host_sbox[pattern[j]];
            hval[0] ^= sbv;
            hval[1] ^= sbv;
        }
        
        uint32_t p;
        if (hashfilter.Get(*((uint64_t*) hval), p) == Ok) {
            assert(p > 0);
        } else
            p = 0;

        if (logbuf_used >= DEF_STRBUF_SIZE)  {
            perror("Not enough memory for logbuf, please make DEF_STRBUF_SIZE larger ");
            exit(-1);
        }

        char* cstr = logbuf + logbuf_used;
        memcpy(cstr, &p, sizeof(uint32_t));
        strcpy(cstr + sizeof(uint32_t), pattern.c_str());

        if (hashfilter.Put(*((uint64_t*) hval), logbuf_used) != Ok) {
            perror("Error while buiding hash table");
            exit(-1);
        }
        logbuf_used += pattern.size() + sizeof(uint32_t);// + 1;
        logbuf_used = (logbuf_used | (MEM_ALIGNMENT - 1)) + 1;
        //printf("logbuf_used = %u\n", logbuf_used);
    }
    inPhrases.close();
   
    hashfilter.BuildBF();

    strbuf = new char[DEF_STRBUF_SIZE];
    strbuf_used = MEM_ALIGNMENT;
    if (!strbuf) {
        perror("cannot allocate host memory for strbuf\n");
        exit(-1);
    }
    memset(strbuf, 0, DEF_STRBUF_SIZE);

    int total = 0;
    for (int i = 0; i < hashfilter.num_buckets; i++) {
        for (int j = 0; j < hashfilter.bucket_size; j++) {
            uint32_t p = hashfilter.ReadValue(i, j);
            if (!p) continue;
            //printf("(%u,%u) %x\n", i, j, p);
            assert(p < logbuf_used);
            assert((p & (MEM_ALIGNMENT - 1)) == 0);


            uint32_t val = strbuf_used;
            strbuf_used += sizeof(uint32_t);
            uint32_t num = 0;
            while (p) {
                num ++;
                char* q = logbuf + p + sizeof(uint32_t);
                int len = strlen(q);
                strcpy(strbuf + strbuf_used, q);
                strbuf_used += len + 1;
                p = *(uint32_t*) ( logbuf + p);

            }
            memcpy(strbuf + val, &num, sizeof(uint32_t));
            strbuf_used = (strbuf_used | (MEM_ALIGNMENT - 1)) + 1;

            hashfilter.WriteValue(i, j, val);
            total += num;

            // printf("(%u,%u) %u strings, val=%u\n", i, j, num, val);
            // int l = 0;
            // while (num > 0) {
            //     char c = strbuf[val + sizeof(uint32_t) + l];
            //     printf("%c", c);
            //     if (c == '\0') {
            //         printf("\n");
            //         num --;
            //     }
            //     l ++;
            // }
        }
    }

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


#ifdef ENABLE_COUNTING
    exit_on_error("cudaMalloc dev_cnt1",
                  cudaMalloc((void **) &dev_cnt1, 
                             sizeof(unsigned int)));


    exit_on_error("cudaMalloc dev_cnt2",
                  cudaMalloc((void **) &dev_cnt2, 
                             sizeof(unsigned int)));

    exit_on_error("cudaMemset dev_cnt1 = 0",
                  cudaMemset(dev_cnt1, 0, sizeof(unsigned int)));

    exit_on_error("cudaMemset dev_cnt2 = 0",
                  cudaMemset(dev_cnt2, 0, sizeof(unsigned int)));

#endif

    delete [] logbuf;
}


void cleanup() 
{
    cudaFree(dev_strbuf);
    cudaFree(dev_buckets);
    cudaFreeHost(pinnedBuf);

#ifdef ENABLE_COUNTING
    cudaFree(dev_cnt1);
    cudaFree(dev_cnt2);
#endif

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
        blocks_x = 256;
        blocks_y = numthreads / (256 * blocksize);
        threads_1d = blocksize;
    } else if (numthreads > blocksize) {
        blocks_x = numthreads / blocksize;
        blocks_y = 1;
        threads_1d = blocksize;
    } else {
        blocks_x = 1;
        blocks_y = 1;
        threads_1d = numthreads;
    }

    unsigned int threads_used = blocks_y * blocks_x * threads_1d;
    numthreads -= threads_used;
      //printf("pick_dim %d %d %d\n", blocks_y, blocks_x, threads_1d);
    dimGrid = dim3(blocks_x, blocks_y);
    dimBlock = dim3(threads_1d);
    return threads_used;
}


inline __device__ void dev_set_bit(int i, unsigned char *bv) {
    //unsigned int *p = (unsigned int *)bv;
    //unsigned int bitMask = 1 << (i & 31);
    //atomicOr(&p[i >> 5], bitMask);
    unsigned int bitMask = 1 << (i & 7);
    bv[i >> 3] |= bitMask;
}

inline __device__ bool dev_is_bit_set(int i, unsigned char *bv) {
    unsigned int bitMask = 1 << (i & 7);
    return (bv[i >> 3] & bitMask);
}

inline bool is_bit_set(int i, unsigned int *bv) {
    unsigned int word = bv[i >> 5];
    unsigned int bitMask = 1 << (i & 31);
    return (word & bitMask);
}


// inline __device__ bool dev_strcmp(unsigned char* text, 
//                                   unsigned char* patt) 
// {
//     int k = 0;
//     while (patt[k] != '\0') {
//         if (text[k] != patt[k])
//             return false;
//         k ++;
//     }
//     return true;
// }


__device__ bool dev_check_str(int i,
                              uint val,
                              unsigned char* text,
                              unsigned char* dev_strbuf)
{
    
    uint num = *(uint*) (dev_strbuf + val);
    unsigned char *q = dev_strbuf + val + sizeof(uint32_t);
    int l = 0;
    int k = 0;
    bool match = true;
    int watch = 210414200;
    
    while (num) {
        if (q[l] == '\0') {
            //if (i == watch)
            //printf("i = %u, num = %u, %c, %d\n", i , num, q[l], match);

            if (match) {
                //printf("match!!\n");
                return true;
            }
            match = true;
            k = 0;
            num --;
        } else {

            if (q[l] != text[k]) {
                match = false;
            }
            // if (i == watch)
            //     printf("i = %u, comparing %c %c, %d\n", i, q[l], text[k], match);
            k ++;
        }
        l ++;
    }
    return false;
}

__global__ void GrepKernel(unsigned char *d_a,
                           unsigned char *dev_buckets,
                           unsigned char *dev_strbuf,
                           unsigned char *dev_pos_bitmap,
                           unsigned int char_offset,
                           uint32_t strbuf_used,
                           unsigned int* dev_cnt1,
                           unsigned int* dev_cnt2)
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
    
    tag =  hval[0] & TAGMASK;
    tag += (tag == 0); 
    i1 = hval[1] & INDEXMASK;
    i2 = (i1 ^ (tag * 0x5bd1e995)) & INDEXMASK;

    // for test
    //dev_set_bit(i, dev_pos_bitmap);

    // don't know why but this turns out slow ...
    // "BucketType* b1= (BucketType*) (dev_buckets + i1 * sizeof(BucketType));"
    // so I read first 24 bytes from bucket i1 as a ulong3
    // including 8-byte tagbits and 16-byte bfbits

    ulong3 v1 = *(ulong3 *) (dev_buckets + i1 * sizeof(BucketType));
    BucketType* b1 = (BucketType*) &v1;
    uint16_t *tagbits1 = (uint16_t*) (b1->tagbits_);

    for (int j = 0; j < bucket_size; j ++) {
        uint32_t tag1 = tagbits1[j] & TAGMASK;

        if (tag1 == tag) {
            int val_offset = (int) offsetof(BucketType, valbits_) + j * sizeof(uint32_t);
            uint p = *(uint*) (dev_buckets + i1 * sizeof(BucketType) + val_offset);
            bool match = dev_check_str(i, p, d_a + i, dev_strbuf);
            // bool match = true;
            if (match) {
                //printf("=====i=%d, set!\n", i);
                dev_set_bit(i, dev_pos_bitmap);
                //printf("%d\n", dev_is_bit_set(i, dev_pos_bitmap));
                return;
            }
        }
    }

#ifdef ENABLE_COUNTING
    atomicAdd(dev_cnt1, 1);
#endif

#ifdef ENABLE_FILTERING
    unsigned char* bfbits1 = b1->bfbits_;
    if (!(dev_is_bit_set(i2 & BFINDEXMASK, bfbits1) &&
          dev_is_bit_set((i2 / num_bfbits) & BFINDEXMASK, bfbits1))) 
        return;
#endif

#ifdef ENABLE_COUNTING
    atomicAdd(dev_cnt2, 1);
#endif

    // read first 8 bytes from bucket i1
    // including 8-byte tagbits 
    ulong1 v2 = *(ulong1 *) (dev_buckets + i2 * sizeof(BucketType));
    BucketType *b2 = (BucketType*) &v2;
    uint16_t *tagbits2 = (uint16_t*) (b2->tagbits_);

    for (int j = 0; j < bucket_size; j ++) {
        uint32_t tag2 = tagbits2[j] & TAGMASK;
        if (tag2 == tag) {   
            int val_offset = (int) offsetof(BucketType, valbits_) + j * sizeof(uint32_t);
            uint p = *(uint*) (dev_buckets + i2 * sizeof(BucketType) + val_offset);
            bool match = dev_check_str(i, p, d_a + i, dev_strbuf);
            //bool match = true;
            if (match) {
                //printf("=====i=%d, set!\n", i);
                dev_set_bit(i, dev_pos_bitmap);
                return;
            }
        }
    }
    return;
}

void process_corpus(char* corpus_file) 
{
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

        //printf("Executing grep on %d\n", size);

        read(fd, pinnedBuf + offset, size);

        exit_on_error("cudaMemcpyAsync",
                      cudaMemcpyAsync(dev_corpus + offset, pinnedBuf + offset, size, cudaMemcpyHostToDevice, streams[i]));

        unsigned int char_offset = 0;
        while (numthreads > 0) {
            unsigned int tu = pick_dim(dimGrid, dimBlock, numthreads, BLOCK_SIZE);
            printf("Executing GrepKernel (%d,%d,%d) @ %u\n", dimGrid.x, dimGrid.y, dimBlock.x, (offset + char_offset));
            GrepKernel<<<dimGrid, dimBlock, 0, streams[i]>>>(dev_corpus, 
                                                             dev_buckets,
                                                             dev_strbuf,
                                                             dev_pos_bitmap, 
                                                             offset + char_offset,
                                                             strbuf_used,
                                                             dev_cnt1,
                                                             dev_cnt2);


            checkReportCudaStatus("GrepKernel");
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

    printf("GrepKernel done\n");

#ifdef ENABLE_COUNTING

    unsigned int host_cnt1, host_cnt2;
    exit_on_error("cudaMemcpy dev_cnt1 to host_cnt1",
                  cudaMemcpy(&host_cnt1, dev_cnt1,
                             sizeof(unsigned int), cudaMemcpyDeviceToHost));

    exit_on_error("cudaMemcpy dev_cnt2 to host_cnt2",
                  cudaMemcpy(&host_cnt2, dev_cnt2,
                             sizeof(unsigned int), cudaMemcpyDeviceToHost));

    printf("total = %u, cnt1 = %u, cnt2 = %u\n", filesize, host_cnt1, host_cnt2);
#endif

}


void print_positions(char *corpus_file)
{
    printf("print positions\n");

    size_t filesize;
    int f = getfile(corpus_file, &filesize);
    if (f == -1) {
        perror(corpus_file);
        exit(-1);
    }
    char *buf = pinnedBuf;

    filesize = min((unsigned long long)filesize, (unsigned long long)FILE_MAX);

    unsigned int *bv = (unsigned int *)malloc(filesize / 8 + 1);
    if (!bv) {
        perror("cannot allocate cpu memory\n");
        exit(-1);
    }

    exit_on_error("cudaMemcpy corpus results to host",
                  cudaMemcpy(bv, dev_pos_bitmap,
                             filesize / 8 + 1, cudaMemcpyDeviceToHost));

    int file_ints = filesize / 32;
    int prev_end = -1;
    int num_matched = 0;
    printf("filesize = %d, file_ints = %d\n", filesize, file_ints);
    for (int i = 0; i < file_ints; i++) {
        //printf("bv[%d] = %d\n", i, bv[i]);
        if (bv[i]) {
	        for (int j = ffs(bv[i]) - 1; j < 32; j++) {
	            int offset = i*32 + j;
	            if (is_bit_set(offset, bv)) {
	                /* Find end of previous line */
	                if (offset > prev_end && buf[offset] != '\n') {
                        char *sol = ((char*)memrchr(buf, '\n', offset));
                        int start_line;
                        if (sol) {
                            start_line = sol - buf;
                        } else {
                            start_line = 0;
                        }

	                    int end_line;
                        char *eol = (char*)memchr(buf + offset, '\n', filesize - offset + 1);
                        end_line = eol - buf;
                        j += end_line - offset;
	                    if (buf[start_line] == '\n') start_line++;
	                    fwrite(buf + start_line, 1, end_line - start_line, stdout);
	                    fputc('\n', stdout);
	                    prev_end = end_line;
                        num_matched ++;
	                }
	            }
	        }
        }
    }

    close(f);

    printf("\n\n%d lines in corpus matched\n", num_matched);
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
    //cout << sizeof(uint) << " " << sizeof(ulong4) << endl;
    //cout << offsetof(BucketType, valbits_) << endl;
    //return 0;

    timing_stamp("start", false);

    setup(patterns_file);

    timing_stamp("setup", false);

    process_corpus(corpus_file);

    timing_stamp("grep corpus", false);

    print_positions(corpus_file);

    timing_stamp("copying out", false);

    cleanup();

    timing_stamp("cleanup done", true);

    timing_report();

    struct cudaDeviceProp cdp;
    cudaGetDeviceProperties(&cdp, 0);
    printf("\ndeviceOverlap = %d\n", cdp.deviceOverlap);
}

