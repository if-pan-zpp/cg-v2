#pragma once
#include "model/NativeStructure.hpp"
#include <string>
#include <vector>
#include <unordered_map>

namespace CG {
    /* This class represents a single chain in <Model>. In particular,
     * we can use it to generate/change conformations, for example when
     * layout out residues in a (partially) disordered chain in a straight
     * line or a self-avoiding walk. */
    class Chain {
    public:
        std::vector<AminoAcid> residues;
        Real3List positions;
        std::vector<NativeStructure> structured_parts;

        /* Morph the chain into a straight-line formation. */
        void into_line(Real3 start, Real3 direction, Real step);

        /* Morph the chain into a self-avoiding walk.
         * TODO: figure out the randomness in the program. */
        void into_saw();

        /* Apply an affine transform to the chain. */
        void apply(RealAffine3 const& aff);
    };
}