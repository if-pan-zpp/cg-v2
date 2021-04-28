#pragma once
#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/Results.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"
using namespace std;

namespace cg::reference {
    using namespace cg::toolkit;

    /* CPC14.pdf, 3.2 */
    class NativeBondAngle: public Force {
    public:
        const Real k = 30.0*eps/(radian*radian);

        NativeBondAngle(PseudoAtoms const &pseudoAtoms,
                        NativeStructure const &ns);

        void compute(Reals3 &forces) override;
        void dumpResults(Results &results) override;

    private:
        PseudoAtoms const &pseudoAtoms;
        NativeStructure const &ns;

        vector<unsigned char> enabled;
        vector<Real> nativeTheta;
        Real energy;
    };
}
