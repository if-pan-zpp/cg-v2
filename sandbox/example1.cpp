#include "loaders/PDBFile.hpp"
#include "loaders/ParameterFile.hpp"
#include "Simulation.hpp"
#include "forces/local/HarmonicTethers.hpp"
#include "forces/local/NativeBondAngle.hpp"
#include "forces/nonlocal/PauliExclusion.hpp"
#include "forces/nonlocal/NativeContacts.hpp"
#include "integrators/LangevinPredictorCorrector.hpp"
#include "reporters/StateReporter.hpp"
using namespace cg::toolkit;
using namespace cg::reference;

int main() {
    PDBFile pdbFile("data/example1/1ubq.pdb");
    ParameterFile paramFile("data/parametersMDCG.txt");
    pdbFile.fullModel.deriveContactsFromAllAtoms(paramFile.parameters);
    auto model = pdbFile.fullModel.reduce();

    for (int traj = 0; traj < 1; ++traj) {
        Simulation sim(model);
        ModelData const &modelData = sim.modelData;
        
        // Create and set integrator
        LangevinPredictorCorrector lpc(sim.delta, sim.modelData.pseudoAtoms, sim.rng);
        sim.integrator = &lpc;

        // Create verlet list
        vector<pair<int, int>> exclusions;
        Neighborhood::Spec spec {
                .minLocalOffset = 3,
                .include4 = true,
                .exclusions = &exclusions,
                .pairTypes = {{"NORMAL", "NORMAL"}},
                .cutoff = 5.0 * angstrom,
                .pad = 15.0 * angstrom
        };
        Neighborhood const &verletList = sim.topology.createNeighborhood(spec);
        
        // Create and attach forces
        NativeBondAngle nba(modelData.pseudoAtoms, modelData.ns);
        sim.attachForce(&nba);
        
        HarmonicTethers ht(modelData.pseudoAtoms, modelData.ns);
        sim.attachForce(&ht);

        PauliExclusion pe(modelData.pseudoAtoms, sim.topology, verletList); 
        sim.attachForce(&pe);

        NativeContacts nc(modelData.pseudoAtoms, modelData.ns, sim.topology, 5.0 * angstrom);
        sim.attachForce(&nc);
        
        // Create and attach reporters
        StateReporter stateRep(modelData.pseudoAtoms, sim.results, sim.delta);
        sim.attachReporter(&stateRep, 200);

        // Run the simulation
        sim.run(100000);
    }
    return 0;
}
