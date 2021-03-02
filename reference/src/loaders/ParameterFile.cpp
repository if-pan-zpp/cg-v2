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
    for (int i = 0; i < 2; ++i)
        getline(file, line);

    vector<string> pairs_order = {"GG", "GP", "GX", "PG", "PP", "PX", "XG",
                                  "XP", "XX"};

    /* Heurestic bond angle potential parameters. */
    getline(file, line); /* Header */
    for (int i = 0; i < pairs_order.size(); ++i) {
        vector<Real> params(7);
        for (int j = 0; j < params.size(); ++j) {
            file >> params[j];
        }

        bond_angle_params[pairs_order[i]] = params;
    }

    /* Heurestic dihedral angle potential parameters. */
    getline(file, line); /* Header */
    for (int i = 0; i < pairs_order.size(); ++i) {
        vector<Real> params(6);
        for (int j = 0; j < params.size(); ++j) {
            file >> params[j];
        }

        dihedral_angle_params[pairs_order[i]] = params;
    }

    /* Specificities. */
    /* Get the order of amino acids. */
    string order, res_name;
    for (int i = 0; i < NUM_AMINOACIDS; ++i) {
        file >> res_name;
        auto code = (char)ResidueName(res_name);
        order.push_back(code);
        specificities[code] = Specificity();
    }
    /* Get the rest of structure data. */
    for (int i = 0; i < NUM_AMINOACIDS; ++i) {
        file >> specificities[order[i]].polarity;
    }
    for (int i = 0; i < NUM_AMINOACIDS; ++i) {
        file >> specificities[order[i]].coordination_number;
    }
    for (int i = 0; i < NUM_AMINOACIDS; ++i) {
        file >> specificities[order[i]].hydrophobic_coordination_number;
    }
    for (int i = 0; i < NUM_AMINOACIDS; ++i) {
        file >> specificities[order[i]].polar_coordination_number;
    }

    /* Amino acid radii. */
    getline(file, line);
    for (int i = 0; i < NUM_AMINOACIDS; ++i) {
        file >> amino_acid_radii[(char)ResidueName(i)];
    }

    /* Pair distances (and optionally MJ). */
    getline(file, line);
    if (load_mj_matrix) {
        mj_matrix = PairMatrix();
    }

#define AMINO_ACID_PAIRS 210
    for (int i = 0; i < AMINO_ACID_PAIRS; ++i) {
        string res1, res2;
        Real dist;
        file >> res1 >> res2 >> dist;

        string code1((char)ResidueName(res1), 1);
        string code2((char)ResidueName(res2), 1);
        ss_minimal_distances[code1 + code2] = dist;

        if (load_mj_matrix) {
            Real energy;
            file >> energy;
            mj_matrix.value()[code1 + code2] = energy;
        }
    }
}
