#include "data/ModelData.hpp"
#include "utils/Units.hpp"
using namespace cg::reference;
using namespace cg::toolkit;
using namespace std;

ModelData::ModelData(Model const &model) {
    size_t n = 0;
    
    for (auto const &id_and_chain : model.chains) {
        Chain const &chain = id_and_chain.second;
        n += chain.positions.cols();
    }

    pseudoAtoms.n = n;
    pseudoAtoms.pos = Reals3::Zero(3, n);
    pseudoAtoms.vel = Reals3::Zero(3, n);

    // TODO: get pseudoatoms' masses
    pseudoAtoms.mass = Reals::Ones(n) * f77mass;

    pseudoAtoms.type = vector<string>(n, "NORMAL");
    pseudoAtoms.typeRanges["NORMAL"] = {0, n};

    // < -10 means undefined
    ns.tether = -11 * Reals::Ones(n);
    ns.bond = -11 * Reals::Ones(n);
    ns.dihedral = -11 * Reals::Ones(n);

    size_t patom_id = 0;
    for (auto const &id_and_chain : model.chains) {
        Chain const &chain = id_and_chain.second;
        
        size_t length = chain.residues.size();

        assert (chain.positions.cols() == length);
        for (size_t i = 0; i < length; ++i) {
            pseudoAtoms.pos.col(patom_id) = chain.positions.col(i);
            pseudoAtoms.chainId.push_back(id_and_chain.first);

            pseudoAtoms.aminoAcidCode.push_back(aaCodeFromName(string(chain.residues[i])));

            patom_id++;
        }


        for (toolkit::NativeStructure const &tns : chain.structuredParts) {
            for (size_t i = 0; i < tns.tether.cols(); ++i) {
                ns.tether(tns.offset + i) = tns.tether(i);
            }

            for (size_t i = 0; i < tns.bond.cols(); ++i) {
                ns.bond(tns.offset + i) = tns.bond(i);
            }

            for (size_t i = 0; i < tns.dihedral.cols(); ++i) {
                ns.dihedral(tns.offset + i) = tns.dihedral(i);
            }

            for (toolkit::NativeStructure::Contact const &contact : tns.contacts) {
                ns.contacts.push_back(NativeStructure::Contact {
                        .residues = contact.residues,
                        .distance = contact.distance
                });
            }
        }
    }

    // TODO: get the rest of pseudoatoms fields

    // TODO: get the rest of NativeStructure fields
} 
