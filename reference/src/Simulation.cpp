#include "Simulation.hpp"
using namespace cg::reference;
using namespace std;

Simulation::Simulation(toolkit::Model const &model, unsigned seed):
    rng(seed),
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

void Simulation::calcForces() {
    forces = Reals3::Zero(3, modelData.pseudoAtoms.n);
    for (Force *f : forceObjects) {
        f -> compute(forces);
    }
}

void Simulation::dumpResults() {
    results.clear();
    for (Force *f : forceObjects) {
        f -> dumpResults(results);
    }
}

void Simulation::run(int max_steps) {
    assert(integrator);

    // TODO: get real temperature value
    modelData.pseudoAtoms.initMovement(rng, 0.35, delta);

    calcForces();
    integrator -> init(forces);
    
    for (int step = 0; step < max_steps; ++step) {
        // TODO: add a possibility to break here

        topology.update();
        calcForces();
        dumpResults();
        integrator -> step(forces);

        // Go through reporters
        for (const pair<Reporter*, int>& rep_and_period: reporters) {
            if (step % rep_and_period.second == 0) {
                rep_and_period.first -> report(step);
            }
        }

    }
}
