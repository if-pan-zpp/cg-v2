#include "tools/param/LegacyParser.hpp"
#include "utils/Text.hpp"
#include <sstream>
using namespace mdk;
using namespace mdk::param;
using namespace std;

static void fetchDefAngleParams(istream &is, Parameters &data) {
    skipLine(is);
    auto ss = lineStream(is);
    for (auto& coeff: data.defAngleParams)
        ss >> coeff;
}

static void fetchAngleParams(istream &is, Parameters &data) {
    skipLine(is);
    for (auto var: variants()) {
        auto ss = lineStream(is);
        for (auto& coeff: data.angleParams[var])
            ss >> coeff;
    }
}

static void fetchDihedralParams(istream &is, Parameters &data) {
    skipLine(is);
    for (auto var: variants()) {
        auto ss = lineStream(is);
        for (auto& coeff: data.dihedralParams[var])
            ss >> coeff;
    }
}

static void fetchSpecificity(istream &is, Parameters &data) {
    auto ss = lineStream(is);
    vector<AminoAcid> order;
    for (string name; ss >> name; ) {
        auto acid = (AminoAcid)(string)name;
        order.push_back(acid);
    }

    ss = lineStream(is);
    for (auto const& acid: order)
        ss >> (int&)data.specificity[acid].polarization;

    ss = lineStream(is);
    for (auto const& acid: order)
        ss >> data.specificity[acid].coordNum;

    ss = lineStream(is);
    for (auto const& acid: order)
        ss >> data.specificity[acid].hydrophobicCoordNum;

    ss = lineStream(is);
    for (auto const& acid: order)
        ss >> data.specificity[acid].polarCoordNum;
}

static void fetchRadii(istream &is, Parameters &data) {
    skipLine(is);

    auto ss = lineStream(is);
    for (auto const& acid: AminoAcid::allAminoAcids)
        ss >> data.radius[acid];
}

static void fetchPairwiseData(istream &is, Parameters &data) {
    string header;
    getline(is, header);

    bool fetchMJ = (header == "amino acid pair distances and energies");
    if (fetchMJ) {
        data.mjMatrix = Parameters::PerPairData<double>();
    }

    int numPairs = AminoAcid::numAminoAcids * (AminoAcid::numAminoAcids+1) / 2;
    for (int i = 0; i < numPairs; ++i) {
        auto ss = lineStream(is);

        string name1, name2;
        ss >> name1 >> name2;
        AminoAcid acid1(name1), acid2(name2);

        double minDist;
        ss >> minDist;
        data.pairwiseMinDist[{acid1, acid2}] = minDist;
        data.pairwiseMinDist[{acid2, acid1}] = minDist;

        if (fetchMJ) {
            double mjEnergy;
            ss >> mjEnergy;

            data.mjMatrix.value()[{acid1, acid2}] = mjEnergy;
            data.mjMatrix.value()[{acid2, acid1}] = mjEnergy;
        }
    }
}

Parameters LegacyParser::read(istream &is) {
    Parameters data;

    fetchDefAngleParams(is, data);
    fetchAngleParams(is, data);
    fetchDihedralParams(is, data);
    fetchSpecificity(is, data);
    fetchRadii(is, data);
    fetchPairwiseData(is, data);

    return data;
}
