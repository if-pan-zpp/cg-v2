#include <iostream>
#include "loaders/PDBFile.hpp"
#include "loaders/SequenceFile.hpp"
#include "loaders/ParameterFile.hpp"
using namespace CG;
using namespace std;

int main() {
    ParameterFile paramfile("data/parametersMDCG.txt", true);
    PDBFile pdbfile("data/example1/1ubq.pdb");
    SequenceFile seqfile("data/example3/glut.txt");

    pdbfile.full_model.derive_contacts_from_all_atoms(paramfile.parameters);
    auto redux = pdbfile.full_model.reduce();
    return EXIT_SUCCESS;
}
