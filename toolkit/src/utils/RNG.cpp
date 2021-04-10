#include "utils/RNG.hpp"
using namespace cg::toolkit;

#ifdef LEGACY_MODE

RNG::RNG(unsigned seed) {
    idum = -(int) seed;
}

Real RNG::ran2() {
    int k, j;
    if(idum <= 0) {
        idum = max(-idum, 1);
        idum2 = idum;
        for(j = NTAB + 7; j>=0; j--) {
            k = idum / IQ1;
            idum = IA1 * (idum - k * IQ1) - k * IR1;
            if(idum < 0)
                idum += IM1;
            if(j < NTAB)
                iv[j] = idum;
        }
        iy = iv[0];
    }
    k = idum / IQ1;
    idum = IA1 * (idum - k * IQ1) - k* IR1;
    if(idum < 0)
        idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if(idum2 < 0)
        idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = idum;
    if(iy < 1)
        iy += IMM1;
    return min(AM * iy, RNMX);
}

Real RNG::uniform() {
    return ran2();
}

Real RNG::normal() {
    Real r1 = ran2();
    Real r2 = ran2();
    return sqrt(-2.0 * log(r1)) * cos(2 * M_PI * r2);
}

#endif
