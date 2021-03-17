#include "loaders/ParameterFile.hpp"
#include <sstream>
#include <fstream>

using namespace cg::toolkit;
using namespace std;

ParameterFile::ParameterFile(std::istream &file, bool loadMjMatrix) {
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
    Parameters::AngleParams preParams;
    for (Index i = 0; i < pairs_order.size(); ++i) {
        vector<Real> params(7);
        for (Index j = 0; j < params.size(); ++j) {
            file >> params[j];
        }

        preParams[pairs_order[i]] = params;
    }
    for (auto amino_acid1: AminoAcid::allCodes) {
        for (auto amino_acid2: AminoAcid::allCodes) {
            auto& loc = parameters.bondAngleParams[{amino_acid1, amino_acid2}];
            if (amino_acid1 != 'G' && amino_acid1 != 'P')
                amino_acid1 = 'X';
            if (amino_acid2 != 'G' && amino_acid2 != 'P')
                amino_acid2 = 'X';
            loc = preParams[{amino_acid1, amino_acid2}];
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

        preParams[pairs_order[i]] = params;
    }
    for (auto aminoAcid1: AminoAcid::allCodes) {
        for (auto aminoAcid2: AminoAcid::allCodes) {
            auto& loc = parameters.dihedralAngleParams[{aminoAcid1, aminoAcid2}];
            if (aminoAcid1 != 'G' && aminoAcid1 != 'P')
                aminoAcid1 = 'X';
            if (aminoAcid2 != 'G' && aminoAcid2 != 'P')
                aminoAcid2 = 'X';
            loc = preParams[{aminoAcid1, aminoAcid2}];
        }
    }

    /* Specificities. */
    /* Get the order of amino acids. */
    string order, resName;
    for (Index i = 0; i < AminoAcid::numAminoAcids; ++i) {
        file >> resName;
        auto code = (char)AminoAcid(resName);
        order.push_back(code);
        parameters.specificities[code] = Parameters::Specificity();
    }
    /* Get the rest of structure data. Unfortunately we need to write out
     * all four loops. */
    for (Index i = 0; i < AminoAcid::numAminoAcids; ++i) {
        file >> parameters.specificities[order[i]].polarity;
    }
    for (Index i = 0; i < AminoAcid::numAminoAcids; ++i) {
        file >> parameters.specificities[order[i]].coordinationNumber;
    }
    for (Index i = 0; i < AminoAcid::numAminoAcids; ++i) {
        file >> parameters.specificities[order[i]].hydrophobicCoordinationNumber;
    }
    for (Index i = 0; i < AminoAcid::numAminoAcids; ++i) {
        file >> parameters.specificities[order[i]].polarCoordinationNumber;
    }

    /* Amino acid radii. */
    getline(file, line); /* End of last one. */
    getline(file, line); /* Header. */
    for (Index i = 0; i < AminoAcid::numAminoAcids; ++i) {
        file >> parameters.aminoAcidRadii[(char)AminoAcid(i)];
    }

    /* Pair distances (and optionally MJ).
     * TODO: load_mj_matrix must actually "be correct", fix it later. */
    getline(file, line); /* End of last one. */
    getline(file, line); /* Header. */
    if (loadMjMatrix) {
        parameters.mjMatrix = Parameters::PairMatrix();
    }

    auto numAminoAcidPairs = AminoAcid::numAminoAcids * (AminoAcid::numAminoAcids + 1);
    for (Index i = 0; i < numAminoAcidPairs; ++i) {
        string res1, res2;
        Real dist;
        file >> res1 >> res2 >> dist;

        auto code1 = (char)AminoAcid(res1);
        auto code2 = (char)AminoAcid(res2);

        parameters.ssMinimalDistances[{code1, code2}] = dist;
        parameters.ssMinimalDistances[{code2, code1}] = dist;

        if (loadMjMatrix) {
            Real energy;
            file >> energy;
            parameters.mjMatrix[{code1, code2}] = energy;
            parameters.mjMatrix[{code2, code1}] = energy;
        }
    }
}

ParameterFile::ParameterFile(const filesystem::path &path,
    bool loadMjMatrix) {
    auto filestream = ifstream(path);
    *this = ParameterFile(filestream, loadMjMatrix);
}
