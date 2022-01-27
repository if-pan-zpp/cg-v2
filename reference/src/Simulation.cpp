#include "Simulation.hpp"
using namespace cg::reference;
using namespace cg::toolkit;
using namespace std;

Simulation::Simulation(toolkit::Model const &model, toolkit::RNG &rng) :
    rng(rng),
    modelData(model),
    topology(modelData.pseudoAtoms, modelData.ns) {
    
    // TODO: get real temperature value
    modelData.pseudoAtoms.initMovement(rng, 0.35 * epsDivkB, delta);
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

    calcForces();
    dumpResults();
    integrator -> init(forces);
    
    for (const pair<Reporter*, int>& rep_and_period: reporters) {
        rep_and_period.first -> report(0);
    }

    for (int step = 1; step <= max_steps; ++step) {
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
