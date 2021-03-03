#include <iostream>
#include <fstream>
#include "loaders/ParameterFile.hpp"

int main() {
    auto file = std::ifstream("data/parametersMDCG.txt");
    auto paramfile = CG::ParameterFile(file, true);
    return EXIT_SUCCESS;
}
