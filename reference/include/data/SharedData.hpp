#pragma once
#include <map>

using std::vector;
using std::map;

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

    struct LJWallSharedData {
        Reals3 pull_ref_pos;
        vector<unsigned> adiabCoeff;
        map<size_t, bool> connected_to_zwall;
    };

    struct WallSharedData {
        size_t upper_n, lower_n, old_upper_n, old_lower_n; // ip1,ip2,icw(1),icw(2)
        pair<Real3, Real3> wall_forces; // (x,y,z-dforce, x,y,z-uforce)
        Real fresist;
        Real fresistperp;
    };

}
