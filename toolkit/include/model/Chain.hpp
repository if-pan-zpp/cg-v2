#pragma once
#include "model/NativeStructure.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include "utils/Units.hpp"
#include "utils/RNG.hpp"

namespace cg::toolkit {
    /* This class represents a single chain in <Model>. In particular,
     * we can use it to generate/change conformations, for example when
     * layout out residues in a (partially) disordered chain in a straight
     * line or a self-avoiding walk. */
    class Chain {
    public:
        std::vector<AminoAcid> residues;
        Real3List positions;
        /* Native structures over the chain. Note that if they
         * overlap, we should prioritize the structures later in the list. */
        std::vector<NativeStructure> structuredParts;

        /* Morph the chain into a straight-line formation. */
        void intoLine(Real3 start, Real3 direction, Real step);

        /* Morph the chain into a self-avoiding walk.
         * TODO: figure out the randomness in the program. */
        void intoSAW(bool dense, bool use_pbc, Real initial_density, Real cutoff, RNG &rng);

        /* Apply an affine transform to the chain. */
        void apply(RealAffine3 const &aff);

    private:
        Real bond_length = 3.8 * angstrom;
        int SAW_attempts_left = 9000;
        
         // for PBC, wouldn't be needed if this class had reference to Topology
        Real3 cell, cell_inv;

    };
    
}
