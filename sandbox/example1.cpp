#include "loaders/PDBFile.hpp"
#include "loaders/ParameterFile.hpp"
#include "Simulation.hpp"
#include "forces/local/HarmonicTethers.hpp"
#include "forces/local/NativeBondAngle.hpp"
#include "forces/nonlocal/PauliExclusion.hpp"
#include "forces/nonlocal/NativeContacts.hpp"
#include "integrators/LangevinPredictorCorrector.hpp"
#include "reporters/StateReporter.hpp"
#include "reporters/PositionDiffReporter.hpp"
using namespace cg::toolkit;
using namespace cg::reference;

int main() {
    PDBFile pdbFile("data/example1/1ubq.pdb");
    ParameterFile paramFile("data/parametersMDCG.txt");
    pdbFile.fullModel.deriveContactsFromAllAtoms(paramFile.parameters);
    auto model = pdbFile.fullModel.reduce();

    for (int traj = 0; traj < 1; ++traj) {
        RNG rng(448u);
        rng.uniform(); // advance rng state once, because cg.f also does it

        // Initialize positions in the model
        Model m = model;
        for (auto &it : m.chains) {
            it.second.intoSAW(false, false, 0.0, 4.56 * angstrom, rng);
        }


        // Create simulation object
        Simulation sim(m, rng);
        ModelData const &modelData = sim.modelData;
        const PseudoAtoms &p = modelData.pseudoAtoms;


        // Create and set integrator
        LangevinPredictorCorrector lpc(sim.delta, sim.modelData.pseudoAtoms, rng);
        sim.integrator = &lpc;

        
        // Create and attach local forces
        NativeBondAngle nba(modelData.pseudoAtoms, modelData.ns);
        sim.attachForce(&nba);
        
        HarmonicTethers ht(modelData.pseudoAtoms, modelData.ns);
        sim.attachForce(&ht);


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

        // Create and attach global forces
        PauliExclusion pe(modelData.pseudoAtoms, sim.topology, verletList); 
        sim.attachForce(&pe);

        NativeContacts nc(modelData.pseudoAtoms, modelData.ns, sim.topology, 5.0 * angstrom);
        sim.attachForce(&nc);
        

        // Create and attach reporters
        StateReporter stateRep(modelData.pseudoAtoms, sim.results, sim.delta);
        sim.attachReporter(&stateRep, 200);


        // Run the simulation
        sim.run(10000);
    }
    return 0;
}
