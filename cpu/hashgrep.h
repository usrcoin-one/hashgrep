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
*/

#ifndef HASH_GREP_H
#define HASH_GREP_H

#include <string>
#include <queue>
#include <vector>
#include <stdint.h>

#include "../hashes/hash_rot_sbox_pre_2.h"
#include "../hashtables/cuckoohashtable.h"

using namespace std;

#ifndef DEF_PATT_LEN
#define DEF_PATT_LEN 19
#endif



class FilterError
{
    const char* const message;

public:
    FilterError(const char* const msg = 0) : message(msg) {}
    const char* const what()
    {
        return message;
    }
};


class Filter
{
    // rolling hash
    hash_rot_sbox_pre_2<DEF_PATT_LEN> hash;
    
    CuckooHashtable<uint64_t, char*, 16, 20> hashfilter;

    //computes rolling hash functions
    void updateHashes(u_int8_t nextChar);
    
public:
    Filter();
    ~Filter();

    void buildFilter(char* phrasesfilename) throw(FilterError);
    //Filter& operator<<(string& pattern) throw(FilterError);
    /* void loadFilterFromFile(int fd) throw(FilterError); */
    /* void saveFilterToFile(int fd) throw(FilterError); */

    void processCorpus(int fd) throw(FilterError);

    void printStatus() {
        cout << hashfilter.Info() << endl;
    }
};

/* Handy crap for filling in templates */
template <int P>
struct kpowerto
{ enum { value = 2654435769U * kpowerto<P-1>::value };};

template<>
struct kpowerto<0>
{ enum { value = 1 }; };

template <int P>
struct fpowerto
{ enum { value = (0xffffffff & (2166136261U * kpowerto<P-1>::value)) };};

template<>
struct fpowerto<0>
{ enum { value = 1 }; };

#endif
