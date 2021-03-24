#include "reporters/StateReporter.hpp"
#include <iostream>
using namespace cg::reference;
using namespace std;

StateReporter::StateReporter(PseudoAtoms const& _pseudoAtoms) :
    pseudoAtoms(_pseudoAtoms) {
}

void StateReporter::report(int step) {
    Real kin_energy = 0;
    Reals3 const& vel = pseudoAtoms.vel;
    for (size_t i = 0; i < pseudoAtoms.n; ++i) {
        kin_energy += 0.5 * pseudoAtoms.mass[i] * vel.col(i).squaredNorm();
    }
    cout << "Step #" << step << ": kin_energy = " << kin_energy << endl;
}
