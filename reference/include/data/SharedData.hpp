#pragma once
using std::vector;

namespace cg::reference {
    struct CoordNumber {
        uint8_t backbone;
        uint8_t sidechain;
        uint8_t ssHydrophobic;
        uint8_t ssPolar;
    };
    
    struct SharedData {
        vector<CoordNumber> coordNumbers;
        vector<unsigned> neiCount;

        SharedData(size_t n) :
            coordNumbers(n, {0, 0, 0, 0}),
            neiCount(n, 0) {
        }
    };
}
