#include "Simulation.hpp"
using namespace cg::reference;
using namespace std;

Simulation::Simulation(toolkit::Model const& model):
    modelData(model),
    topology(modelData.pseudoAtoms, modelData.ns) {
    
}

void Simulation::attachForce(Force *force) {
    forceObjects.push_back(force);
}

void Simulation::attachReporter(Reporter *reporter, int period) {
    assert (period > 0);
    reporters.emplace_back(reporter, period);
}

void Simulation::calc_forces() {
    Real energy = 0;

    forces = Reals3::Zero(3, modelData.pseudoAtoms.n);
    for (Force *f : forceObjects) {
        f -> compute(energy, forces);
    }
}

void Simulation::run(int max_steps) {
    assert(integrator);
    
    for (int step = 0; step < max_steps; ++step) {
        // TODO: add a possibility to break here

        calc_forces();
        integrator -> step(delta, forces);

        // Go through reporters
        for (const pair<Reporter*, int>& rep_and_period: reporters) {
            if (step % rep_and_period.second == 0) {
                rep_and_period.first -> report(step);
            }
        }
    }
}
