#pragma once
#include <mdk/utils/AminoAcid.hpp>
#include <mdk/tools/cmap/ContactMap.hpp>
#include <mdk/tools/model/CoarseModel.hpp>

namespace mdk::seq {
    class Sequence {
    public:
        struct Chain {
            std::vector<AminoAcid> codes;
            std::vector<std::string> contactMaps;
        };

        std::optional<double> screeningDist;
        std::vector<Chain> chains;
        std::unordered_map<std::string, cmap::ContactMap> contactMaps;

    public:
        CoarseModel asModel() const;
    };
}