#include "forces/local/HarmonicTethers.hpp"
using namespace cg;

void HarmonicTethers::compute(Reals *energy, Reals3 *force) {
    for (auto [start, end]: ls->chains) {
        for (int ix = start; ix+1 < end; ++ix) {
            auto dist0 = ns->isNative[ix] ? ns->tether(ix) : default_dist0;
            auto tetherVec = system->pos.col(ix+1)-system->pos.col(ix);
            auto tetherVecLen = tetherVec.norm();
            auto diff = tetherVecLen-dist0;
            auto diffsq = diff*diff;

            /* |F| = k_1 x + k_3  x^3; k_1 is H1 and k_3 is H2 from what I see */
            // auto en = k1*diffsq/2.0 + k3*diffsq*diffsq/4.0;
            auto fvalue = k1*diff + k3*diff*diffsq;

            /* As for the force cap, why should it even be here? */

            tetherVec /= tetherVecLen;
            force->col(ix) -= tetherVec * fvalue;
            force->col(ix+1) += tetherVec * fvalue;
        }
    }
}
