#include "loaders/PDBFile.hpp"
#include "loaders/ParameterFile.hpp"
#include "Simulation.hpp"
#include "forces/local/HarmonicTethers.hpp"
#include "forces/local/LocalRepulsive.hpp"
#include "forces/local/NativeBondAngle.hpp"
#include "forces/local/ComplexNativeDihedral.hpp"
#include "forces/nonlocal/NativeContacts.hpp"
#include "forces/nonlocal/PauliExclusion.hpp"
#include "integrators/LangevinPredictorCorrector.hpp"
#include "reporters/PositionDiffReporter.hpp"
#include "reporters/StateReporter.hpp"
#include "data/NativeStructure.hpp"
using namespace cg::toolkit;
using namespace cg::reference;

int main() {
    PDBFile pdbFile("data/1ubq.pdb");
    ParameterFile paramFile("data/parametersMDCG.txt");
    pdbFile.fullModel.deriveContactsFromAllAtoms(paramFile.parameters);
    auto model = pdbFile.fullModel.reduce();

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
    PseudoAtoms const &p = modelData.pseudoAtoms;
    
    // Create and set integrator
    LangevinPredictorCorrector lpc(sim.delta, sim.modelData.pseudoAtoms, rng);
    sim.integrator = &lpc;

    // Create and attach forces
    HarmonicTethers ht(modelData.pseudoAtoms, modelData.ns);
    sim.attachForce(&ht);

    LocalRepulsive lr(modelData.pseudoAtoms);
    sim.attachForce(&lr);

    NativeBondAngle nba(modelData.pseudoAtoms, modelData.ns);
    sim.attachForce(&nba);

    ComplexNativeDihedral cnd(modelData.pseudoAtoms, modelData.ns);
    sim.attachForce(&cnd);

    NativeContacts nc(modelData.pseudoAtoms,
                      modelData.ns,
                      sim.topology,
                      18.0 * angstrom);
    sim.attachForce(&nc);

    vector<pair<int, int>> nativeContacts;
    for (cg::reference::NativeStructure::Contact const &cnt : modelData.ns.contacts) {
        nativeContacts.emplace_back(cnt.residues);
    }
    sort(nativeContacts.begin(), nativeContacts.end());
    
    Neighborhood::Spec spec {
        .minLocalOffset = 3,
        .include4 = true,
        .exclusions = &nativeContacts,
        .pairTypes = {{"NORMAL", "NORMAL"}},
        .cutoff = 5.0 * angstrom,
        .pad = 15.0 * angstrom
    };
    Neighborhood const &verletList = sim.topology.createNeighborhood(spec);

    PauliExclusion pe(modelData.pseudoAtoms,
                      sim.topology,
                      verletList);
    sim.attachForce(&pe);

    // Create and attach reporter
    PositionDiffReporter diffRep(modelData.pseudoAtoms, "data/positions.txt");
    sim.attachReporter(&diffRep, 10);

    // StateReporter sr(modelData.pseudoAtoms, sim.results, sim.delta);
    // sim.attachReporter(&sr, 10);

    // Run the simulation
    sim.run(500);
    return 0;
}
