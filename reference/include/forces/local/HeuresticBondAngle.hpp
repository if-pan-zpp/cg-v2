#pragma once
#include <unordered_map>
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"
#include "utils/TupleHash.hpp"

namespace cg::reference {
    /* This is the bond angle potential for unstructured parts,
     * see 3.2.2 in CPC14.pdf */
    class HeuresticBondAngle: public Force {
    public:
        enum ResidueType { GLY, PRO, X };
        using Params = std::unordered_map<std::pair<ResidueType, ResidueType>, Real>;

    private:
        PseudoAtoms const *pseudoAtoms;
        // Needed for checking whether part is unstructured.
        NativeStructure const *ns;
        Params params;

    public:
        HeuresticBondAngle(PseudoAtoms const& pseudoAtoms, NativeStructure const& ns, Params const& params);

        void compute(Real &energy, Reals3 &forces) override;
    };
}