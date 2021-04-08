#pragma once

#include "forces/Force.hpp"
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include "utils/Units.hpp"
#include "data/Results.hpp"
#include <vector>
using namespace std;

namespace cg::reference {
    using namespace cg::toolkit;

    /* CPC14.pdf, 3.2.1, different for simple and complex
     * variants because of the hyperparameters. */
    class ComplexNativeDihedral : public Force {
    public:
        Real K1 = 0.33 * eps/(radian*radian);
        Real K3 = 0.33 * eps/(radian*radian);

        ComplexNativeDihedral(PseudoAtoms const &pseudoAtoms, NativeStructure const &ns);

        void compute(Reals3 &forces) override;
        void dumpResults(Results &results) override;
    private:
        PseudoAtoms const &pseudoAtoms;

        Real energy;
        vector<unsigned char> enabled;
        vector<Real> nativePhi;
    };
}
