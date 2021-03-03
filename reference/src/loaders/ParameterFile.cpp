#include "loaders/ParameterFile.hpp"
#include "utils/ResidueName.hpp"
#include <sstream>

using namespace CG;
using namespace std;

ParameterFile::ParameterFile(std::istream &file, bool load_mj_matrix) {
    string line;

    /* First two lines are about "simple bending angle potential parameters",
     * but they initialize sigma1 (cg.f:5523), which in turn is all over the
     * place; we shall therefore for now ignore it.
     * TODO: figure out what is this. */
    for (Index i = 0; i < 2; ++i)
        getline(file, line);

    vector<string> pairs_order = {"GG", "GP", "GX", "PG", "PP", "PX", "XG",
                                  "XP", "XX"};

    /* Heurestic bond angle potential parameters. */
    getline(file, line); /* Header */
    for (Index i = 0; i < pairs_order.size(); ++i) {
        vector<Real> params(7);
        for (Index j = 0; j < params.size(); ++j) {
            file >> params[j];
        }

        bond_angle_params[pairs_order[i]] = params;
    }

    /* Heurestic dihedral angle potential parameters. */
    getline(file, line); /* Rest of last line */
    getline(file, line); /* Header */
    for (Index i = 0; i < pairs_order.size(); ++i) {
        vector<Real> params(6);
        for (Index j = 0; j < params.size(); ++j) {
            file >> params[j];
        }

        dihedral_angle_params[pairs_order[i]] = params;
    }

    /* Specificities. */
    /* Get the order of amino acids. */
    string order, res_name;
    for (Index i = 0; i < NUM_AMINOACIDS; ++i) {
        file >> res_name;
        auto code = (char)ResidueName(res_name);
        order.push_back(code);
        specificities[code] = Specificity();
    }
    /* Get the rest of structure data. */
    for (Index i = 0; i < NUM_AMINOACIDS; ++i) {
        file >> specificities[order[i]].polarity;
    }
    for (Index i = 0; i < NUM_AMINOACIDS; ++i) {
        file >> specificities[order[i]].coordination_number;
    }
    for (Index i = 0; i < NUM_AMINOACIDS; ++i) {
        file >> specificities[order[i]].hydrophobic_coordination_number;
    }
    for (Index i = 0; i < NUM_AMINOACIDS; ++i) {
        file >> specificities[order[i]].polar_coordination_number;
    }

    /* Amino acid radii. */
    getline(file, line); /* End of last one. */
    getline(file, line); /* Header. */
    for (Index i = 0; i < NUM_AMINOACIDS; ++i) {
        file >> amino_acid_radii[(char)ResidueName(i)];
    }

    /* Pair distances (and optionally MJ). */
    getline(file, line); /* End of last one. */
    getline(file, line); /* Header. */
    if (load_mj_matrix) {
        mj_matrix = PairMatrix();
    }

#define AMINO_ACID_PAIRS 210
    for (Index i = 0; i < AMINO_ACID_PAIRS; ++i) {
        string res1, res2;
        Real dist;
        file >> res1 >> res2 >> dist;

        string code1(1, (char)ResidueName(res1));
        string code2(1, (char)ResidueName(res2));
        auto pair1 = code1 + code2;
        auto pair2 = code2 + code1;
        ss_minimal_distances[pair1] = dist;
        ss_minimal_distances[pair2] = dist;

        if (load_mj_matrix) {
            Real energy;
            file >> energy;
            mj_matrix.value()[pair1] = energy;
            mj_matrix.value()[pair2] = energy;
        }
    }
}
