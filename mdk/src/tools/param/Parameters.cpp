#include "tools/param/Parameters.hpp"
using namespace mdk::param;
using namespace std;

std::vector<Variants> mdk::param::variants() {
    return { GG, GP, GX, PG, PP, PX, XG, XP, XX };
}

Parameters::Parameters() {
    defAngleParams = Coeffs(7);
    for (auto const& var: variants()) {
        angleParams[var] = defAngleParams;
        dihedralParams[var] = defAngleParams;
    }

    for (auto const& acid1: AminoAcid::allAminoAcids) {
        specificity[acid1] = SpecificityParams();
        radius[acid1] = 0;

        for (auto const& acid2: AminoAcid::allAminoAcids) {
            pairwiseMinDist[{acid1, acid2}] = 0;
        }
    }
}
