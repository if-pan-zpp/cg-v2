#include "loaders/ParameterFile.hpp"
#include "utils/AminoAcid.hpp"
#include <sstream>
#include <fstream>

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

    /* In the parameter file, amino acids other than G and P are not
     * distinguished. */
    vector<Parameters::AminoAcidPair> pairs_order = {
            {'G', 'G'}, {'G', 'P'}, {'G', 'X'}, {'P', 'G'}, {'P', 'P'},
            {'P', 'X'}, {'X', 'G'}, {'X', 'P'}, {'X', 'X'}
    };

    /* Heurestic bond angle potential parameters. */
    getline(file, line); /* Header */

    /* In final bond_angle_params we store all pairs rather than just the equiv
     * classes. Thus we need to proceed in a two-staged fashion. */
    Parameters::AngleParams pre_params;
    for (Index i = 0; i < pairs_order.size(); ++i) {
        vector<Real> params(7);
        for (Index j = 0; j < params.size(); ++j) {
            file >> params[j];
        }

        pre_params[pairs_order[i]] = params;
    }
    for (auto amino_acid1: AminoAcid::all_codes) {
        for (auto amino_acid2: AminoAcid::all_codes) {
            auto& loc = parameters.bond_angle_params[{amino_acid1, amino_acid2}];
            if (amino_acid1 != 'G' && amino_acid1 != 'P')
                amino_acid1 = 'X';
            if (amino_acid2 != 'G' && amino_acid2 != 'P')
                amino_acid2 = 'X';
            loc = pre_params[{amino_acid1, amino_acid2}];
        }
    }

    /* Heurestic dihedral angle potential parameters. */
    getline(file, line); /* Rest of last line */
    getline(file, line); /* Header */
    for (Index i = 0; i < pairs_order.size(); ++i) {
        vector<Real> params(7);
        for (Index j = 0; j < params.size(); ++j) {
            file >> params[j];
        }

        pre_params[pairs_order[i]] = params;
    }
    for (auto amino_acid1: AminoAcid::all_codes) {
        for (auto amino_acid2: AminoAcid::all_codes) {
            auto& loc = parameters.dihedral_angle_params[{amino_acid1, amino_acid2}];
            if (amino_acid1 != 'G' && amino_acid1 != 'P')
                amino_acid1 = 'X';
            if (amino_acid2 != 'G' && amino_acid2 != 'P')
                amino_acid2 = 'X';
            loc = pre_params[{amino_acid1, amino_acid2}];
        }
    }

    /* Specificities. */
    /* Get the order of amino acids. */
    string order, res_name;
    for (Index i = 0; i < AminoAcid::num_amino_acids; ++i) {
        file >> res_name;
        auto code = (char)AminoAcid(res_name);
        order.push_back(code);
        parameters.specificities[code] = Parameters::Specificity();
    }
    /* Get the rest of structure data. Unfortunately we need to write out
     * all four loops. */
    for (Index i = 0; i < AminoAcid::num_amino_acids; ++i) {
        file >> parameters.specificities[order[i]].polarity;
    }
    for (Index i = 0; i < AminoAcid::num_amino_acids; ++i) {
        file >> parameters.specificities[order[i]].coordination_number;
    }
    for (Index i = 0; i < AminoAcid::num_amino_acids; ++i) {
        file >> parameters.specificities[order[i]].hydrophobic_coordination_number;
    }
    for (Index i = 0; i < AminoAcid::num_amino_acids; ++i) {
        file >> parameters.specificities[order[i]].polar_coordination_number;
    }

    /* Amino acid radii. */
    getline(file, line); /* End of last one. */
    getline(file, line); /* Header. */
    for (Index i = 0; i < AminoAcid::num_amino_acids; ++i) {
        file >> parameters.amino_acid_radii[(char)AminoAcid(i)];
    }

    /* Pair distances (and optionally MJ).
     * TODO: load_mj_matrix must actually "be correct", fix it later. */
    getline(file, line); /* End of last one. */
    getline(file, line); /* Header. */
    if (load_mj_matrix) {
        parameters.mj_matrix = Parameters::PairMatrix();
    }

    auto num_amino_acid_pairs = AminoAcid::num_amino_acids * (AminoAcid::num_amino_acids + 1);
    for (Index i = 0; i < num_amino_acid_pairs; ++i) {
        string res1, res2;
        Real dist;
        file >> res1 >> res2 >> dist;

        auto code1 = (char)AminoAcid(res1);
        auto code2 = (char)AminoAcid(res2);

        parameters.ss_minimal_distances[{code1, code2}] = dist;
        parameters.ss_minimal_distances[{code2, code1}] = dist;

        if (load_mj_matrix) {
            Real energy;
            file >> energy;
            parameters.mj_matrix[{code1, code2}] = energy;
            parameters.mj_matrix[{code2, code1}] = energy;
        }
    }
}

ParameterFile::ParameterFile(const filesystem::path &path,
    bool load_mj_matrix) {
    auto filestream = ifstream(path);
    *this = ParameterFile(filestream, load_mj_matrix);
}
