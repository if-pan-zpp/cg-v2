#include <mdk/tools/seq/LegacyParser.hpp>
#include <mdk/tools/pdb/Parser.hpp>
#include <fstream>
using namespace mdk;
using namespace std;

int main() {
    auto model = seq::LegacyParser().read("data/9AAC.txt").asModel();

    auto rand = FortranRandom(0);
    auto density = 1e-4*atom/pow(angstrom, 3);
    auto ixDist = 4.56*angstrom;
    model.morphIntoSAW(rand, false, density, ixDist);

    ofstream coarseFile("data/9AAC.pdb");
    pdb::Parser().write(coarseFile, pdb::Data(AtomicModel(model)));

    return 0;
}