#include "model/Model.hpp"
using namespace CG;
using namespace std;

Model::Model(const Chain &chain) {
    chains[0] = chain;
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

    return *this;
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

void Model::derive_contacts_from_CA_atoms(const Parameters &parameters) {
    for (auto const& [i1, chain1]: chains) {
        auto nresidues1 = chain1.residues.size();
        for (auto const& [i2, chain2]: chains) {
            auto nresidues2 = chain2.residues.size();
            for (Index j1 = 0; j1 < nresidues1; ++j1) {
                auto& v1 = chain1.positions.col(j1);
                Index j2 = (i1 == i2 ? j1 + 3 : 0);
                for (; j2 < nresidues2; ++j2) {
                    auto& v2 = chain2.positions.col(j2);
                    auto dist = (v2 - v1).norm();
                    if (dist < parameters.native_contact_cutoff) {
                        contacts.push_back((Contact) {
                            .res1 = {i1, j1},
                            .res2 = {i2, j2},
                            .distance = dist,
                            .type = "bb" /* Not sure about how to call it. */
                        });
                    }
                }
            }
        }
    }
}
