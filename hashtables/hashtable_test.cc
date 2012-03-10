
#include <random>

#include "cuckoohashtable.h"


using namespace std;
std::mt19937 engine;

int
main() {

    const size_t idxbits = 10;

    CuckooHashtable<uint64_t, char*, 16, idxbits> HT;

    for (size_t i = 0; i  < 4 * (1<< idxbits); i ++) {
        uint64_t key;
        uint32_t* kbuf= (uint32_t*) &key;
        kbuf[0] = engine();
        kbuf[1] = engine();

        //printf("data = 0x%x %x\n", kbuf[0], kbuf[1]);
        char* val = new char[i+1];
        printf("val = 0x%p\n", val);
        if (HT.Put(key, val) != Ok) {
            cout << i << " can not put" << endl;
            break;
        }
        char* newval;
        if (HT.Get(key, newval) != Ok) {
            cout << "can not get" << endl;
            break;
        }
        printf("newval = 0x%p\n", newval);
    }
    cout << HT.Info() << endl;

}
