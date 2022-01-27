#pragma once

#include "utils/Types.hpp"
#include <random>
using namespace std;

#ifndef LEGACY_MODE
namespace cg::toolkit {
    class RNG {
    public:
        RNG(unsigned seed) :
            engine(seed),
            uniform_dist(0., 1.),
            normal_dist(0., 1.) {}

        inline Real uniform() { // ~U(0, 1)
            return uniform_dist(engine);
        }
        inline Real normal() {  // ~N(0, 1)
            return normal_dist(engine);
        }

    private:
        mt19937 engine;
        uniform_real_distribution<Real> uniform_dist;
        normal_distribution<Real> normal_dist;
    };
}

#else // Legacy version used for testing
namespace cg::toolkit {
    class RNG {
    public:
        RNG(unsigned seed);
        Real uniform(); // ~U(0, 1)
        Real normal();  // ~N(0, 1)

    private:
        Real ran2();

        static const int
            IM1 = 2147483563,
            IM2 = 2147483399,
            IMM1 = IM1-1,
            IA1 = 40014,
            IA2 = 40692,
            IQ1 = 53668,
            IQ2 = 52774,
            IR1 = 12211,
            IR2 = 3791,
            NTAB = 32,
            NDIV = 1+IMM1/NTAB;

        static constexpr Real
            EPS = 1.2e-7,
            RNMX = 1. - EPS,
            AM = (float) 1./IM1;

        int iy = 0,
            idum2 = 123456789,
            idum = -448;
        int iv[32];
    };
}
#endif
