#include "reporters/PositionDiffReporter.hpp"
#include "utils/Units.hpp"
#include <iostream>
#include <iomanip>
using namespace std;
using namespace cg::reference;

PositionDiffReporter::PositionDiffReporter(PseudoAtoms const &pseudoAtoms,
                                           string filePath):
    pseudoAtoms(pseudoAtoms) {

    ifstream input(filePath, ifstream::in);

    string token;
    unsigned step_nr;

    while (input >> token) {
        if (token != "STEP") continue;
        input >> step_nr;

        Reals3 pos = Reals3::Zero(3, pseudoAtoms.n);
        for (size_t i = 0; i < pseudoAtoms.n; ++i) {
            input >> pos.col(i)[0] >> pos.col(i)[1] >> pos.col(i)[2];
        }

        referencePositions[step_nr] = pos;
    }
}

void PositionDiffReporter::report(int step) {
    auto it = referencePositions.find(step);
    if (it != referencePositions.end()) {
        Reals3 diff = (it -> second) - (pseudoAtoms.pos / toolkit::angstrom);
        auto diff_len = diff.colwise().norm();
        Real max_diff = diff_len.maxCoeff();
        Real mean_diff = diff_len.mean();

        cout << "Step #" << setw(8) << step
             << ": max " << setw(12) << max_diff 
             << ": mean " << setw(12) << mean_diff << endl;
    }
    else {
        cerr << "No information about positions in step " << step << endl;
    }
}
