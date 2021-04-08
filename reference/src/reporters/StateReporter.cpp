#include "reporters/StateReporter.hpp"
#include <iostream>
#include <iomanip>
using namespace cg::reference;
using namespace std;

StateReporter::StateReporter(PseudoAtoms const& _pseudoAtoms,
                             Results const &_results,
                             Real _delta) :
    pseudoAtoms(_pseudoAtoms),
    results(_results),
    delta(_delta) {
}

void StateReporter::report(int step) {
    Real kinEnergy = 0;
    Reals3 const& vel = pseudoAtoms.vel;
    for (size_t i = 0; i < pseudoAtoms.n; ++i) {
        kinEnergy += 0.5 * pseudoAtoms.mass[i] * vel.col(i).squaredNorm();
    }
    kinEnergy /= delta * delta;

    Real totEnergy = kinEnergy + results.potEnergy;
    cout << "Step #" << setw(6) << step
         << ": kinEnergy = " << setw(8) << kinEnergy
         << ": totEnergy = " << setw(8) << totEnergy
         << endl;
}
