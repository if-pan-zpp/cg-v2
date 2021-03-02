#include <iostream>
#include <fstream>
#include "geometry/PDBFile.hpp"

int main() {
    auto pdbfile = std::ifstream("1ubq.pdb");
    auto pdb = geometry::PDBFile(pdbfile);
    return EXIT_SUCCESS;
}
