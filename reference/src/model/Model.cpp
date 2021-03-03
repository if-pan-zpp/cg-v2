#include "model/Model.hpp"
using namespace CG;
using namespace std;

Model::Model(const Chain &chain, const string &name) {
    chains[name] = chain;
}

Model &Model::operator+=(const Model &model2) {
    /* Determine safe offset distance for new chains. */
    Index offset = 0;
    for (auto const& [ix, chain]: chains) {
        offset = max(offset, ix);
    }

    /* Insert the chains. */
    for (auto const& [ix, chain]: model2.chains) {
        chains[offset + ix] = chain;
    }
}

Model Model::operator+(const Model &model2) const {
    auto sum = *this;
    sum += model2;
    return sum;
}

void Model::apply(RealAffine3 const& aff) {
    for (auto& [name, chain]: chains) {
        chain.apply(aff);
    }
}
