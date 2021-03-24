#include "reporters/StateReporter.hpp"
#include <iostream>
using namespace cg::reference;
using namespace std;

StateReporter::StateReporter(PseudoAtoms const& _pseudoAtoms) :
    pseudoAtoms(_pseudoAtoms) {
}

void StateReporter::report(int step) {
    Real kin_energy = 0;
    // TODO: calculate kinetic energy
    cout << "Step #" << step << ": kin_energy = " << kin_energy << endl;
}
