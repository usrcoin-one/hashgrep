/* -*- Mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*- */
#ifndef _HASH_ROT_SBOX_PRE_2_H_
#define _HASH_ROT_SBOX_PRE_2_H_

#include "hash_common.h"
#include "hash_buf.h"

template<int LEN> class hash_rot_sbox_pre_2 {
private:
    static const uint64_t sbox[];
    uint64_t b[LEN];
    int bufptr;



    bool full;

    static const uint32_t shift3 = (3 * LEN % 32);

public:

    hash_rot_sbox_pre_2()
    {
    }

    uint64_t h;

    inline uint32_t hval() { return (uint32_t)(h >> 32); }
    inline uint32_t hval1() { return (uint32_t)h; }
    inline uint32_t hval2() { return (uint32_t)(h >> 32); }

    inline void reset() {
        h = 0;
        bufptr = 0;
        for (int i = 0; i < LEN; i++) {
            b[i] = 0;
        }
        full = false;
    }
    inline void init() {
        reset();
    }

    inline void update(unsigned char c) {
        h = rol64(h, 1);
        h ^= sbox[c];
        h ^= b[bufptr];
        b[bufptr++] = rol64(sbox[c], LEN);
        if (bufptr >= LEN) {
            bufptr = 0;
            full = true;
        }
    }

    inline uint8_t is_full()
    {
        return full;
    }
};

template <int LEN> const uint64_t hash_rot_sbox_pre_2<LEN>::sbox[] = {
     0x198c9c06f9ed325b, 0x98fdeeaca8935b2b, 0xb54edd4194c5926e, 0x4b3e5aff2251c980,
     0x46a2b7aeda0083cd, 0x1955a0ce10a87a3d, 0xc2d772897315f773, 0x6b234d9215415bc1,
     0x67cfc73e49a0155c, 0x76b46e3b5d1411b0, 0x8dc9a383c9b2bb51, 0xbaa5d992bc9aff2f,
     0xf53d3b4cad2663cf, 0xa58a6ac26221adbe, 0x2bfa5b92ef9f5d5, 0x2124661ec274a5d5,
     0x969bcbeae5ddb10, 0x562a1f042bcd908f, 0x359f604d1be85524, 0xbdfdefeee13b67a6,
     0x6de664d4ee6387b9, 0x50635ebb850883b1, 0x3371158f3a9224b9, 0x2588ae5f33dc8a6e,
     0xac34275c3596ca5b, 0x27fee973168be4b6, 0xc6ce09dcab3fc56b, 0x8a901564776fcfaa,
     0x159e0019c036c75e, 0x248a696dd6fa76a1, 0xc21fae6380ae9d69, 0xa5b26570eeb567d5,
     0xc5adaf96f5db7761, 0xac5e3ec00e75aec8, 0x3a9de0a7f875a979, 0xd67f1ce2d89513a7,
     0xac35381824a0e6fa, 0xd8add854b414fddb, 0x66f4920d2b402af2, 0xd708caf803a2b300,
     0x5c386396b9e788b3, 0xcd2b393ed715c898, 0x3bd7bed308a21e4a, 0x7c2ca3651a85ddbd,
     0x45f7ad5c3f4fa6d5, 0x91920b16cfcb10fa, 0x2a7dbf333f2bc06d, 0x811760b43a64478e,
     0x8a5b2a4526e6e8c1, 0x1b189b7417cbda38, 0xcb346b8e394cca62, 0x5cedfba5da675169,
     0x4e5fd67568b7be42, 0x999382b31b4253fc, 0xf51d8e9b1a845de4, 0xd4402fc84aa7bb01,
     0xa4db67c95020023c, 0x3228b24e2da5917c, 0xb62489c6cf54a57a, 0x6133cd88e9eda7a4,
     0xc0f0935b341f5a5e, 0x9bc4d42f158d6456, 0x1fed6b06a3cebc43, 0x615ce06f434e24a6,
     0xdc43539d85d8ae30, 0xf512b3aa97b6bea5, 0x15d4e4c7f5373368, 0xde78cadc705aa659,
     0x98c9c7cfa094a6eb, 0xb69edabb369f37b9, 0x3b4436984016c913, 0x519fbbc2c11904cb,
     0x31affb8bdce4a035, 0x399aad9aa67f8ce2, 0x70e6ae021fe0fcd0, 0xf288b2cf389be9f4,
     0xf9ce3401d98914ce, 0xa9bb8e4410d90af8, 0x59011343b3f4ba5b, 0xd43607ac7d876750,
     0xc10f2c7cba19c68d, 0x958cc99f558be83e, 0x63897faa42384ac, 0x57a678d386255d27,
     0x45a6d6fca3d6b39b, 0xa880f682aa719902, 0xc4b648e294699b51, 0xfa8b567ae6dffd53,
     0x20dba6635e13b77a, 0x2b52fafb12ef1934, 0x2494ea2c50af07bd, 0x216ab3a9fa8d7902,
     0x4077e9c681d89fe4, 0xae4c0eb762d87b95, 0xabe58dea777ea66a, 0x11f00622c6a1afe,
     0xe65ab596fdf73318, 0x4779229a867b5085, 0xcc343c7090f25f3e, 0x574800ce2cd2601d,
     0x2ac4f76ab6aace23, 0x14bc41c7a0d2ab05, 0x6bdc1678071a63be, 0x5d6a8fe4b0f9db4b,
     0xdfbb095fed7a92f3, 0xdea5de5574509bec, 0x8252c8f44ea9104b, 0x58089c88bb37dc0d,
     0x54b49501dc0f6607, 0x75ce75fc2f70abdf, 0xd2fe32885b7a141b, 0x7d87f655a3bd9269,
     0x945b728eb22bf7ee, 0xb8504da2c1085122, 0x1e585d3508573c06, 0x9fa8eac04c0b7a93,
     0x4debde69eed5bb11, 0x250acedf0f7ebe06, 0x34f3e6c4043ead38, 0xdeba5719bdd960a,
     0x36186dee5e2f8c51, 0x1cbc4978844109e1, 0x22aff2de810c5682, 0x81906afc240b4b78,
     0xb65ec442980e4e68, 0x22a7145d15acf230, 0x6a7eaf355292f41f, 0xbc2bac9523c22972,
     0x62ba4d390361aff2, 0xc71ebf0e6ef668e1, 0x946902ea625e94f, 0x9368110461399186,
     0x70ad260a82437675, 0xc05332c066d2e225, 0xf2e859e1c3232dcc, 0xf02ecceaa6d66d79,
     0x6dc4557a8d6ddc39, 0x814544bf63ac17db, 0xbdaeee77193c4113, 0x7cb47d631b77b8bb,
     0xdbaf9cb55c7ea9b, 0xdbb3d142f0c0ffb1, 0x2b0e6fac9de103b7, 0x3b4f4403a88fa5e9,
     0xa25b044d85ebd9ec, 0x4a8e523f3c4b0eec, 0x488a6d045955edc7, 0x2c32eea06d0e4994,
     0xe51c101a53e28045, 0x6fcb00dd24d2fcf3, 0xad0161c49aec1855, 0x97679a0bd48f0bf2,
     0x34c780ac05a5621a, 0xbce0227c67959463, 0xd128b891c5c2a19e, 0x23b1d6dfb33c146,
     0x81c5308878bc064f, 0x2603cf1319a6c555, 0x5019d3fd801ac7ff, 0x3aeb7155c6dc626a,
     0xd35611fb512d99c7, 0x7e5cb80ee834ddec, 0x8ed635bb15c1830c, 0x514f201d2e770bfb,
     0x64e8e1886c1f119b, 0xbcced7bea676e276, 0xfe0f5b2e2f1ac4ca, 0xac27fa35aa26aff8,
     0xa027055c8bc33c39, 0x9dc60274c199b353, 0x2f34120eff663c1f, 0xc232821359493df1,
     0x4f114121c159329c, 0x4fa2902b764a08a8, 0x5fcbb325d37fb2d, 0x14395f0b158aa327,
     0x5cfa050b223c3c2a, 0xc780c7799ff9df67, 0xc7730bac85006a, 0x3c159a9e24b283f0,
     0x6bd833608e070044, 0x8bd7f37608ce44b3, 0x77d8c52fff59b439, 0xbf3fbfbba9cf5fcf,
     0x637c2a3d14253ff7, 0x7a8ae3b6d52b9b18, 0xea93e0166c2c42c6, 0x78017e79ea6dcf66,
     0xda86cec46d9b24e4, 0x8dd09fb5e8fe5c43, 0x6cf6f13a13a030f7, 0xd886cc042f2da725,
     0xe9629b20d41beaf9, 0x775931aa1f286d55, 0xae350711e6538976, 0xa84edf8975cd70b6,
     0xabc9a739945af361, 0x1ab7a53ce695c7de, 0x17edc22e1767c02b, 0xd4ed70209c285fad,
     0xf894aa7fdf27b31, 0xdddf6b3d6d9d043e, 0x1183dc794dcab59a, 0x967359ba71f92b3,
     0x1b13a7cd3db68d00, 0x9556718b1198519e, 0x8ad633c0a2bc44bc, 0x94be88df32f0219d,
     0x96d168da6c424dfb, 0x16274c142585b712, 0x28b748cd8c09db3c, 0xa405949c37412b63,
     0x7e6298c57fd83820, 0xd8ffcdae4b99d012, 0x7a2c11d31880c032, 0xf87fec913ca22e54,
     0x37ea65ccf688b58b, 0x61f96333962ef8c1, 0xb89ce510ba718cb2, 0xfd4fcbc654315caf,
     0x6be12f5f6ef45e62, 0xbe33688ef899370f, 0xd965c821389bebf6, 0x66593593d4e42d71,
     0x429d0c6d62d9cf73, 0x215de6388f46b125, 0xbe7d6d9503f369b, 0xa092a95de1d523f9,
     0x63b5b7834fe670ce, 0xabf2c8b89b9cd1be, 0x1c7ea180d64a8e60, 0x676f0d9bf51f2fb4,
     0xc58bd4ac7f7616b3, 0x573bb39f27a9e34, 0x152e54f670269dfa, 0xbf81e21b407f3f74,
     0x4101c6785b6a64bf, 0x7a11a0b78652c96c, 0x5c36ee8d07369ae, 0x722967f08dcafef0,
     0xcb4b018ee32050b2, 0x9d9eb062ec4efd70, 0xcdfc02984e9fda59, 0xbe6035d8bbeae47e,
     0x797514fa1be020b5, 0xe31be6267ab0fcb1, 0x64e017862acdf5af, 0xa2849b94601aa2e6,
     0xc1430b519190e30a, 0xf0e67ef18ed275e3, 0x4f808a0b1abb56a3, 0x3aacd69b4c871c0f,
     0x5a27af8ab1a45d33, 0xe03ce6dc90e5fcd4, 0xcb2c318d54ebb53c, 0x6d943f9cbe2423c6
};

#endif /* _HASH_ROT_SBOX_PRE_2_H_ */
