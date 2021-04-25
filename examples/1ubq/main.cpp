#include <mdk/tools/pdb/Parser.hpp>
#include <fstream>
#include <iostream>
using namespace mdk;
using namespace std;

int main() {
    ifstream pdb_file("data/1ubq.pdb");
    auto model = pdb::Parser().read(pdb_file).asModel();
    model.addNativeContacts();

    auto coarsened = model.coarsen();

    auto rand = FortranRandom(0);
    auto density = 1e-4*atom/pow(angstrom, 3);
    auto ixDist = 4.56*angstrom;
    coarsened.morphIntoSAW(rand, false, density, ixDist);

    ofstream coarse_file("data/1ubq.saw.pdb");
    pdb::Parser().write(coarse_file, pdb::Data(AtomicModel(coarsened)));

    return 0;
}