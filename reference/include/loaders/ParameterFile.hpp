#pragma once
#include <istream>
#include <vector>
#include <optional>
#include "math/Types.hpp"

namespace CG {
    class ParameterFile {
    public:
        explicit ParameterFile(std::istream& file, bool load_mj_matrix = false);

    private:
        /* We use an alphabet with G (Gly), P (Pro) and X (rest) */
        using AngleParams = std::unordered_map<std::string, std::vector<Real>>;
        AngleParams bond_angle_params;
        AngleParams dihedral_angle_params;

        /* Used in QA potential. See README.txt for more info. */
        struct Specificity {
            /* 0 for Gly/Pro, 1 for hydrophobic, 2 for polar,
             * 4 and 5 for charged. */
            int polarity;

            /* These limit the total number of contacts in QA potential.
             * See CPC14.pdf. */
            int coordination_number;
            int hydrophobic_coordination_number;
            int polar_coordination_number;
        };
        std::unordered_map<char, Specificity> specificities;

        /* Amino acid radii. */
        std::unordered_map<char, Real> amino_acid_radii;

        /* ss-type contact minimal distances for all pairs. */
        using PairMatrix = std::unordered_map<std::string, Real>;
        PairMatrix ss_minimal_distances;

        /* (Optionally) MJ matrix. */
        std::optional<PairMatrix> mj_matrix;
    };
}