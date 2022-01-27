#include "reporters/StateReporter.hpp"
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <iomanip>
using namespace cg::reference;
using namespace cg::toolkit;
using namespace std;

StateReporter::StateReporter(PseudoAtoms const &_pseudoAtoms,
                             Results const &_results,
                             Real _delta) :
    pseudoAtoms(_pseudoAtoms),
    results(_results),
    delta(_delta) {
}

void StateReporter::report(int step) {
    Real kinEnergy = 0;
    Reals3 const &vel = pseudoAtoms.vel;
    for (size_t i = 0; i < pseudoAtoms.n; ++i) {
        kinEnergy += 0.5 * pseudoAtoms.mass[i] * vel.col(i).squaredNorm();
    }
    kinEnergy /= delta * delta;

    Real totEnergy = kinEnergy + results.potEnergy;
    cout << "Step #" << setw(6) << step
         << ": kinEnergy = " << setw(8) << kinEnergy
         << ": totEnergy = " << setw(8) << totEnergy
         << ": actContacts = " << results.activeContacts
         << ": asphericity = " << calculateAsphericity()
         << endl;
}

double StateReporter::calculateAsphericity(void) {
    Reals3 const &pos = pseudoAtoms.pos;
    Real3 middle;
    for (size_t i = 0; i < pseudoAtoms.n; ++i) {
        middle += pos.col(i);
    }
    middle /= pseudoAtoms.n;
    
    Real3_3 mit;
    Real crg;
    for (size_t i = 0; i < pseudoAtoms.n; ++i) {
        Real3 diff = pos.col(i) - middle;
        crg += diff.squaredNorm();
        Real3 sq_diff = diff.array() * diff.array();
        mit(0,0) += sq_diff(1) + sq_diff(2);
        mit(1,1) += sq_diff(2) + sq_diff(0);
        mit(2,2) += sq_diff(0) + sq_diff(1);
        mit(0,1) -= diff(0) * diff(1);
        mit(1,2) -= diff(1) * diff(2);
        mit(0,2) -= diff(0) * diff(2);
    }
    mit(1,0) = mit(0,1);
    mit(2,1) = mit(1,2);
    mit(2,0) = mit(0,2);
    crg = sqrt(crg / pseudoAtoms.n);
    
    auto eigenvaluesC = mit.eigenvalues();
    Real eigenvalues[3] = {
      sqrt(eigenvaluesC(0).real() / pseudoAtoms.n),
      sqrt(eigenvaluesC(1).real() / pseudoAtoms.n),
      sqrt(eigenvaluesC(2).real() / pseudoAtoms.n)
    };
    
    sort(eigenvalues, eigenvalues + 3);
    return eigenvalues[1] - 0.5 * (eigenvalues[0] + eigenvalues[2]);
}
