#include "model/Chain.hpp"

using namespace cg::toolkit;

void Chain::intoLine(Real3 start, Real3 direction, Real step) {
    direction = direction.normalized() * step;
    for (Index i = 0; i < positions.cols(); ++i) {
        positions.col(i) = start + direction * i;
    }
}

void Chain::intoSAW() {
    /* TODO: implement */
}

void Chain::apply(RealAffine3 const &aff) {
    for (Index i = 0; i < positions.cols(); ++i) {
        positions.col(i) = aff * positions.col(i);
    }
}
