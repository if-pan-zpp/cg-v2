#pragma once
#include <vector>
#include <unordered_map>
#include <optional>
#include <mdk/utils/AminoAcid.hpp>
#include <mdk/utils/TupleHash.hpp>

namespace mdk::param {
    enum Variants { GG, GP, GX, PG, PP, PX, XG, XP, XX };
    std::vector<Variants> variants();

    class Parameters {
    public:
        Parameters();

        using Coeffs = std::vector<double>;
        Coeffs defAngleParams;

        template<typename T>
        using PerVariantData = std::unordered_map<Variants, T>;

        PerVariantData<Coeffs> angleParams;
        PerVariantData<Coeffs> dihedralParams;

        struct SpecificityParams {
            enum { GLY_PRO, HYDROPHOBIC, POLAR, POS, NEG } polarization;
            int coordNum, hydrophobicCoordNum, polarCoordNum;
        };

        template<typename T>
        using PerAcidData = std::unordered_map<AminoAcid, T>;

        PerAcidData<SpecificityParams> specificity;
        PerAcidData<double> radius;

        template<typename T>
        using PerPairData = std::unordered_map<
            std::pair<AminoAcid, AminoAcid>, T>;

        PerPairData<double> pairwiseMinDist;
        std::optional<PerPairData<double>> mjMatrix;
    };
}