#include "utils/Random.hpp"
using namespace mdk;
using namespace std;

double Random::uniform(double a, double b) {
    return a + (b - a) * uniform();
}

double Random::normal(double mu, double sigma) {
    return mu + sigma * normal();
}

double ModernRandom::uniform() {
    return unifDist(eng);
}

double ModernRandom::normal() {
    return normalDist(eng);
}

FortranRandom::FortranRandom(int seed) {
    idum = -seed;
    for (auto& x: iv) x = 0;
}

double FortranRandom::uniform() {
    int k, j;
    if (idum <= 0) {
        idum2 = idum = max(-idum, 1);
        for (j = ntab + 7; j >= 0; --j) {
            k = idum / iq1;
            idum = ia1 * (idum - k * iq1) - k * ir1;
            if (idum < 0) idum += im1;
            if (j < ntab) iv[j] = idum;
        }
    }

    k = idum / iq1;
    idum = ia1 * (idum - k * iq1) - k * ir1;
    if (idum < 0) idum += im1;

    k = idum2 / iq2;
    idum2 = ia2 * (idum2 - k * iq2) - k * ir2;
    if (idum2 < 0) idum2 += im2;

    j = iy / ndiv;
    iy = iv[j] - idum2;
    iv[j] = idum;
    if (iy < 1) iy += imm1;

    return min(am * iy, rnmx);
}
