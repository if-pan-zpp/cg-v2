#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"
using namespace std;

namespace cg::reference {
    using namespace cg::toolkit;

    /* CPC14.pdf, 3.2 */
    class NativeBondAngle: public Force {
    public:
        Real k = 30.0*eps/(radian*radian);

        NativeBondAngle(PseudoAtoms const& pseudo_atoms, NativeStructure const& ns);

        void compute(Real &energy, Reals3 &forces) override;

    private:
        PseudoAtoms const& pseudo_atoms;
        NativeStructure const& ns;

        vector<unsigned char> enabled;
        vector<Real> native_theta;
    };
}
