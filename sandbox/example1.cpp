#include "loaders/PDBFile.hpp"
#include "loaders/ParameterFile.hpp"
#include "Simulation.hpp"
#include "forces/local/HarmonicTethers.hpp"
#include "forces/local/NativeBondAngle.hpp"
#include "integrators/LangevinPredictorCorrector.hpp"
#include "reporters/StateReporter.hpp"
using namespace cg::toolkit;
using namespace cg::reference;

int main() {
    PDBFile pdbFile("data/example1/1ubq.pdb");
    ParameterFile paramFile("data/parametersMDCG.txt");
    pdbFile.fullModel.deriveContactsFromAllAtoms(paramFile.parameters);
    auto model = pdbFile.fullModel.reduce();

    for (int traj = 0; traj < 3; ++traj) {
        Simulation sim(model);
        ModelData const& modelData = sim.modelData;
        

        // Create and set integrator
        LangevinPredictorCorrector lpc(sim.modelData.pseudoAtoms);
        sim.integrator = &lpc;


        // Create and attach forces
        NativeBondAngle nba(modelData.pseudoAtoms, modelData.ns);
        sim.attachForce(&nba);

        HarmonicTethers ht(modelData.pseudoAtoms, modelData.ns);
        sim.attachForce(&ht);


        // Create and attach reporters
        StateReporter stateRep(modelData.pseudoAtoms);
        sim.attachReporter(&stateRep, 200);


        // Run the simulation
        sim.run(10000);
    }
    return 0;
}
