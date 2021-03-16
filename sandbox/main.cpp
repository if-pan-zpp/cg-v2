#include <iostream>
#include "loaders/PDBFile.hpp"
#include "loaders/SequenceFile.hpp"
#include "loaders/ParameterFile.hpp"
using namespace cg;
using namespace std;

int main() {
    ParameterFile paramfile("data/parametersMDCG.txt", true);
    PDBFile pdbfile("data/example1/1ubq.pdb");
    SequenceFile seqfile("data/example3/glut.txt");

    pdbfile.fullModel.deriveContactsFromAllAtoms(paramfile.parameters);
    auto redux = pdbfile.fullModel.reduce();
    return EXIT_SUCCESS;
}
