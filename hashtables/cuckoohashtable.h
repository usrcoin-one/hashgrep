/* -*- Mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*- */
#ifndef _CUCKOO_HASHTABLE_H_
#define _CUCKOO_HASHTABLE_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include <sstream>
#include <string>

#include <stdint.h>

using namespace std;

enum Status {
    Ok = 0,
    NotFound = 1,
    NotEnoughSpace = 2,
    NotSupported = 3,
};

template <typename KeyType, typename ValueType, size_t num_tagbits, size_t num_indexbits>
class CuckooHashtable {

    static const size_t bucket_size = 4;

    static const size_t num_keybits = sizeof(KeyType) << 3;
    static const size_t num_valbits = sizeof(ValueType) << 3;
    static const size_t num_buckets = 1ULL << num_indexbits;

    static const size_t INDEXMASK = num_buckets - 1;
    static const uint32_t TAGMASK = (1ULL << num_tagbits) - 1;
    static const uint32_t VALMASK = (1ULL << num_valbits) - 1;
     
    struct Bucket {
        char tagbits_[num_tagbits * bucket_size / 8];
        ValueType valbits_[bucket_size];
    }  __attribute__((__packed__));

    Bucket* buckets_;
    size_t num_keys_;
    static const size_t max_num_keys_ = bucket_size * num_buckets;

    static const size_t MAX_CUCKOO_COUNT = 128;

    inline size_t IndexHash(const char* pkey) const {
        assert(sizeof(KeyType) >= 8);
        uint32_t *data= (uint32_t*) pkey;
        uint32_t r = data[1] & INDEXMASK;
        return  r;
    }

    inline uint32_t TagHash(const char* pkey) const {
        assert(sizeof(KeyType) >= 8);
        uint32_t *data= (uint32_t*) pkey;
        uint32_t r =  data[0] & TAGMASK;
        return r + (r == 0);         // make sure taghash will never be 0
    }

    inline size_t AltIndex(const size_t index, const uint32_t tag) const {
        // 0x5bd1e995 is the hash constant from MurmurHash2
        return (index ^ (tag * 0x5bd1e995)) & INDEXMASK;
    }

    inline uint32_t ReadTag(const size_t i, const size_t j) {
        size_t offset = num_tagbits * j;
        uint32_t v = *(uint32_t *) (buckets_[i].tagbits_ + offset / 8);
        v = (v >> (offset & 0x7)) & TAGMASK;
        return v;
    }

    inline void WriteTag(const size_t i, const size_t j, uint32_t tag) {
        size_t offset = num_tagbits * j;
        uint32_t* p = (uint32_t *) (buckets_[i].tagbits_ + offset / 8);
        *p |= (tag & TAGMASK) << (offset & 0x7);
    }

    inline ValueType ReadValue(const size_t i, const size_t j) const {
        return buckets_[i].valbits_[j];
    }

    inline void WriteValue(const size_t i, const size_t j, ValueType val) {
        buckets_[i].valbits_[j] = val;
    }


public:

    explicit CuckooHashtable(): num_keys_(0) { 
        buckets_ = new Bucket[num_buckets];
    }

    ~CuckooHashtable() { 
        delete [] buckets_;
    }

    // Add a key to the filter.
    Status Put(const KeyType key, const ValueType val);

    // Report if the key is inserted, with false positive rate.
    Status Get(const KeyType key, ValueType& val) const;

    /* methods for providing stats  */
    // summary infomation
    string Info() const {
        stringstream ss;
        ss << "CuckooHashtable Status:" << endl;
        ss << "\t\tNumber of index bits: " << num_indexbits << endl;
        ss << "\t\tNumber of tag bits: " << num_tagbits << endl;
        ss << "\t\tKeys stored: " << Size() << endl;
        ss << "\t\tLoad facotr: " << LoadFactor() << endl;
        ss << "\t\tHashtable size: " << (SizeInBytes() >> 10) << " KB\n";
        if (Size() > 0)
            ss << "\t\tbit/key:   " << BitsPerKey() << endl;
        else
            ss << "\t\tbit/key:   N/A" << endl;

        return ss.str();
    }

    // number of current inserted keys;
    size_t Size() const {return num_keys_;}

    // size of the filter in bytes.
    size_t SizeInBytes() const { return sizeof(Bucket) * num_buckets; }

    double LoadFactor() const { return 1.0 * num_keys_  / max_num_keys_; }

    double BitsPerKey() const { return 8.0 * SizeInBytes() / num_keys_; }

}; // declaration of class CuckooHashtable



template <typename KeyType, typename ValueType, size_t num_tagbits, size_t num_indexbits>
Status 
CuckooHashtable<KeyType, ValueType, num_tagbits, num_indexbits>::Put(const KeyType key, const ValueType v) {

    size_t i = IndexHash((char*) &key);
    uint32_t tag = TagHash((char*) &key);
    ValueType val = v;
   
    for (uint32_t count = 0; count < MAX_CUCKOO_COUNT; count ++) {
        //printf("cout = %u, i = %zu \n", count, i);
        for (size_t j = 0; j < bucket_size; j ++) {

            if (ReadTag(i, j) == 0) {
                WriteTag(i, j, tag);
                WriteValue(i, j, val);
                num_keys_ ++;
                return Ok;
            }
        }

        if (count > 0) {
            size_t j = rand() % bucket_size;
            uint32_t oldtag = ReadTag(i, j);
            ValueType oldval = ReadValue(i, j);
            WriteTag(i, j, tag);
            WriteValue(i, j, val);

            tag = oldtag;
            val = oldval;

        }

        i = AltIndex(i, tag);
    }

    return NotEnoughSpace;
}


template <typename KeyType, typename ValueType, size_t num_tagbits, size_t num_indexbits>
Status 
CuckooHashtable<KeyType, ValueType, num_tagbits, num_indexbits>::Get(const KeyType key, ValueType& val) const {
    bool found = false;
    size_t i1, i2;

    uint32_t tag = TagHash((char*) &key);
    i1 = IndexHash((char*) &key);
    i2 = AltIndex(i1, tag);

    for (size_t j = 0; j < bucket_size; j ++) {

        if (ReadTag(i1, j) == tag) {
            val = ReadValue(i1, j);
            return Ok;
        }

        if (ReadTag(i2, j) == tag) {
            val = ReadValue(i2, j);
            return Ok;
        }        
    }

    return NotFound;
}


#endif // #ifndef _CUCKOO_HASHING_H_
