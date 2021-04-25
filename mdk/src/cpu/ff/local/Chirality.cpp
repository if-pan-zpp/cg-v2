#include "cpu/ff/local/Chirality.hpp"
using namespace mdk;

Chirality::Chirality(const Particles &pts): LocalFF<Chirality>(pts) {
    this->pts = &pts;
    Cnat = Scalars(pts.n);
    d0CubeInv = Scalars(pts.n);

    for (auto const& [start, end]: pts.chains) {
        for (int i = start+1; i < end-2; ++i) {
            auto r0 = vec(pts.pos, i-1), r1 = vec(pts.pos, i),
                 r2 = vec(pts.pos, i+1), r3 = vec(pts.pos, i+2);
            auto v0 = r1 - r0, v1 = r2 - r1, v2 = r3 - r2;

            d0CubeInv[i] = pow(v0.norm(), -3.0);
            Cnat[i] = v0.cross(v1).dot(v2) * d0CubeInv[i];
        }
    }
}

void Chirality::localTerm(double &V, Vectors &gradV, int i, int beg, int end) {
    if (i <= beg || i > end-2) return;

    auto r0 = vec(pts->pos, i-1), r1 = vec(pts->pos, i),
        r2 = vec(pts->pos, i+1), r3 = vec(pts->pos, i+2);
    auto v0 = r1 - r0, v1 = r2 - r1, v2 = r3 - r2;

    auto Ci = v0.cross(v1).dot(v2) * d0CubeInv[i];
    V += 0.5 * echi * pow(Ci - Cnat[i], 2.0);

    auto s = echi * (Ci - Cnat[i]) * d0CubeInv[i];
    auto c01 = s*v0.cross(v1), c02 = s*v0.cross(v2), c12 = s*v1.cross(v2);

    vec(gradV, i-1) = -c12;
    vec(gradV, i) = c12+c02;
    vec(gradV, i+1) = -c02-c01;
    vec(gradV, i+2) = c01;
}
