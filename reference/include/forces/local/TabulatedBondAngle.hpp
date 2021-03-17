#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"

namespace cg::reference {
    /* CPC14.pdf, 3.2.3 */
    class TabulatedBondAngle: public Force {
    public:
        using Params = Eigen::Matrix<Real, 6, Eigen::Dynamic>;

    private:
        PseudoAtoms const *pseudoAtoms;
        NativeStructure const *ns;
        Params params;

    public:
        TabulatedBondAngle(PseudoAtoms const& pseudoAtoms, NativeStructure const& ns, Params params);

        void compute(Real& energy, Reals3& forces) override;
    };
}