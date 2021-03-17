#include "loaders/PDBFile.hpp"
#include "loaders/ParameterFile.hpp"
#include "data/ModelData.hpp"
#include "forces/local/HarmonicTethers.hpp"
#include "forces/local/SimpleNativeDihedral.hpp"
#include "forces/local/NativeBond.hpp"
#include "forces/nonlocal/NativeContacts.hpp"
#include "forces/nonlocal/PauliExclusion.hpp"
#include "integrators/LangevinPredictorCorrector.hpp"
#include "reporters/StateReporter.hpp"
using namespace cg::toolkit;
using namespace cg::reference;

int main() {
    PDBFile pdbFile("data/example1/1ubq.pdb");
    ParameterFile paramFile("data/parametersMDCG.txt");
    pdbFile.fullModel.deriveContactsFromAllAtoms(paramFile.parameters);
    auto model = pdbFile.fullModel.reduce();

    ModelData md(model);
    Topology top(md.pseudoAtoms, md.ns);
    LangevinPredictorCorrector lpc;
    StateReporter sr(md.pseudoAtoms);

    SimpleNativeDihedral snd(md.pseudoAtoms, md.ns);
    lpc.attachForce(&snd);

    NativeBond nb(md.pseudoAtoms, md.ns);
    lpc.attachForce(&nb);

    HarmonicTethers ht(md.pseudoAtoms, md.ns);
    lpc.attachForce(&ht);

    NativeContacts nc(md.ns, top);
    lpc.attachForce(&nc);

    PauliExclusion pe(md.pseudoAtoms, top);
    lpc.attachForce(&pe);

    for (int traj = 0; traj < 10; ++traj) {
        lpc.step(15000);
        sr.report();
    }

    return 0;
}