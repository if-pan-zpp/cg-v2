#include "forces/local/NativeBond.hpp"
using namespace cg;

void NativeBond::compute(Reals *energy, Reals3 *force) {
    for (auto [start, end]: ns->nativeParts) {
        for (int ix = start+1; ix+1 < end; ++ix) {
            auto v0 = system->pos.col(ix) - system->pos.col(ix-1);
            auto v0_norm = v0.norm();

            auto v1 = system->pos.col(ix+1) - system->pos.col(ix);
            auto v1_norm = v1.norm();

            auto theta = acos(-v0.dot(v1)/(v0_norm*v1_norm));

            /* CPC14.pdf has a type on the potential term for the bond
             * angle. */
            // auto en = k*(theta-ns->bond(ix))**2;
            auto fvalue = 0.0;

            /* Here forces would be added */
        }
    }
}