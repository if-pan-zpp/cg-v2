#include "data/ModelData.hpp"
using namespace cg::reference;
using namespace cg::toolkit;

ModelData::ModelData(Model const& model) {
    size_t n = 0;
    
    for (auto const& id_and_chain : model.chains) {
        Chain const& chain = id_and_chain.second;
        n += chain.positions.cols();
    }

    pseudoAtoms.n = n;
    pseudoAtoms.pos = Reals3::Zero(3, n);
    pseudoAtoms.vel = Reals3::Zero(3, n);

    // TODO: get pseudoatoms' masses
    pseudoAtoms.mass = Reals::Ones(n);

    ns.tether = -Reals::Ones(n); // -1 means undefined tether

    size_t patom_id = 0;
    for (auto const & id_and_chain : model.chains) {
        Chain const& chain = id_and_chain.second;
        
        for (size_t i = 0; i < chain.positions.cols(); ++i) {
            pseudoAtoms.pos.col(patom_id) = chain.positions.col(i);
            ns.chainId.push_back(id_and_chain.first);

            patom_id++;
        }


        for (toolkit::NativeStructure const& tns : chain.structuredParts) {
            for (size_t i = 0; i < tns.tether.cols(); ++i) {
                ns.tether(tns.offset + i) = tns.tether(i);
            }
        }
    }

    // TODO: get the rest of pseudoatoms fields

    // TODO: get the rest of NativeStructure fields
} 