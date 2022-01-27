#include "model/FullModel.hpp"
#include "utils/Angles.hpp"
#include <unordered_set>
using namespace cg::toolkit;
using namespace std;

FullModel &FullModel::operator+=(const FullModel &fullModel2) {
    auto offset = chains.size();

    /* Load chains. */
    for (auto const &[chainId, chain]: fullModel2.chains) {
        chains[offset + chainId] = chain;
    }

    /* Load contacts, need to fix offsets. */
    for (auto contact: fullModel2.contacts) {
        contact.res1.first += offset;
        contact.res2.first += offset;
        contacts.push_back(move(contact));
    }

    return *this;
}

FullModel FullModel::operator+(const FullModel &fullModel2) const {
    auto sum = *this;
    sum += fullModel2;
    return sum;
}

void FullModel::apply(const RealAffine3 &aff) {
    for (auto &[chainId, chain]: chains) {
        for (auto &residue: chain) {
            for (auto &[name, pos]: residue.atoms) {
                pos = aff * pos;
            }
        }
    }
}

Model FullModel::reduce() {
    Model redux;

    /* Reduce the chains. */
    for (auto const &[chainId, chain]: chains) {
        redux.chains[chainId] = reduceChain(chain, intraChainContacts[chainId]);
    }

    /* Contacts stay essentially the same. */
    for (auto const &contact: contacts) {
        redux.contacts.push_back((Model::Contact) {
            .res1 = contact.res1,
            .res2 = contact.res2,
            .distance = contact.distance,
            .type = contact.type
        });
    }

    return redux;
}

cg::toolkit::Chain FullModel::reduceChain(
    FullModel::Chain const &chain,
    unordered_set<IndexPair> const &intraContacts) const {

    cg::toolkit::Chain reduxChain;
    NativeStructure ns;

    ns.offset = 0;

    Real3List CA;
    resizeVectorList(CA, chain.size());

    /* Load CA atom positions. */
    for (Index i = 0; i < chain.size(); ++i) {
        CA.col(i) = chain.at(i).atoms.at("CA");
        reduxChain.residues.push_back(chain[i].type);
    }

    /* Derive tether distances. */
    ns.tether.resize(chain.size());
    for (Index i = 0; i+1 < chain.size(); ++i) {
        ns.tether(i) = (CA.col(i+1) - CA.col(i)).norm();
    }

    /* Derive bond angles. */
    ns.bond.resize(chain.size());
    for (Index i = 1; i+1 < chain.size(); ++i) {
        ns.bond(i) = bond(CA.col(i-1), CA.col(i), CA.col(i+1));
    }

    /* Derive dihedral angles. */
    ns.dihedral.resize(chain.size());
    for (Index i = 2; i+1 < chain.size(); ++i) {
        ns.dihedral(i) = dihedral(CA.col(i-2), CA.col(i-1), CA.col(i), CA.col(i+1));
    }

    /* Store intra chain contacts. */
    for (auto const &[i, j] : intraContacts) {
        ns.contacts.emplace_back(NativeStructure::Contact {
                .residues = make_pair(i, j),
                .distance = (CA.col(j) - CA.col(i)).norm()
        });
    }

    reduxChain.positions = move(CA);
    reduxChain.structuredParts.push_back(ns);
    return reduxChain;
}

void FullModel::deriveContactsFromAllAtoms(const Parameters &parameters) {
    std::unordered_set<std::string> backbone = {"N", "CA", "C", "O"};
    Real alpha = pow(26.0/7.0, 1.0/6.0); /* cg.f:4906 */

    /* As for the monstrosity below, unfortunately C++ doesn't have generators
     * (yet, anyways). */

    /* Over chains. */
    for (auto const &[i1, chain1]: chains) {
        for (auto const &[i2, chain2]: chains) {
            if (i2 > i1) continue;

            /* Over residues. */
            for (Index j1 = 0; j1 < chain1.size(); ++j1) {
                auto &res1 = chain1[j1];

                Index j2 = (i1 == i2 ? j1 + 3 : 0);
                for (; j2 < chain2.size(); ++j2) {
                    auto &res2 = chain2[j2];

                    /* Over atoms. */
                    for (auto const &[name1, pos1]: res1.atoms) {
                        if (parameters.atomRadii.count({res1.type, name1}) == 0)
                            continue;
                        auto radius1 = parameters.atomRadii.at({res1.type, name1});
                        char type1 = backbone.find(name1) == backbone.end() ? 'b' : 's';

                        for (auto const &[name2, pos2]: res2.atoms) {
                            if (parameters.atomRadii.count({res2.type, name2}) == 0)
                                continue;
                            auto radius2 = parameters.atomRadii.at({res2.type, name2});
                            char type2 = backbone.find(name2) == backbone.end() ? 'b' : 's';

                            /* Check if in contact. */
                            auto dist = (pos2 - pos1).norm();
                            if (dist < alpha * (radius1 + radius2)) {
                                contacts.push_back((Contact) {
                                    .res1 = {i1, j1},
                                    .res2 = {i2, j2},
                                    .distance = dist,
                                    .type = string(1, type1) + string(1, type2),
                                });

                                if (i1 == i2) {
                                    intraChainContacts[i1].insert({j1, j2});
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
