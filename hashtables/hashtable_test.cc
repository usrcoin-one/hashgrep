
#include <random>

#include "cuckoohashtable.h"


using namespace std;
std::mt19937 engine;

int
main() {

    const size_t idxbits = 20;
    auto HT = new CuckooHashtable<uint64_t, uint32_t, 16, idxbits>;

    for (size_t i = 0; i  < 4 * (1<< idxbits); i ++) {
        uint64_t key;
        uint32_t* kbuf= (uint32_t*) &key;
        kbuf[0] = engine();
        kbuf[1] = engine();

        //printf("data = 0x%x %x\n", kbuf[0], kbuf[1]);
        uint32_t val = i;
        if (HT->Put(key, val) != Ok) {
            cout << "can not put" << endl;
            break;
        }
    }
    cout << HT->Info() << endl;

}
