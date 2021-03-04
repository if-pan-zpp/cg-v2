#pragma once
#include <vector>
#include "utils/Types.hpp"
#include "utils/AminoAcid.hpp"
#include "utils/TupleHash.hpp"

namespace CG {
    /* This class contains extra parameters (like amino acid radii), usually
     * loaded from parameters.txt.
     * TODO: decide whether all these should be in a single data structure. */
    class Parameters {
    public:
        Parameters();

        /* Parameters for heurestic angle force fields. */
        using AminoAcidPair = std::pair<AminoAcid, AminoAcid>;
        using AngleParams = std::unordered_map<AminoAcidPair, std::vector<Real>>;
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
        std::unordered_map<AminoAcid, Specificity> specificities;

        std::unordered_map<AminoAcid, Real> amino_acid_radii;

        /* ss-type contact minimal distances for all pairs. */
        using PairMatrix = std::unordered_map<AminoAcidPair, Real>;
        PairMatrix ss_minimal_distances;

        PairMatrix mj_matrix;

        /* Atom radii.
         * Note: we don't set it in parameters.txt, though perhaps it should
         * be settable externally. */
        using AtomID = std::pair<AminoAcid, std::string>;
        std::unordered_map<AtomID, Real> atom_radii;

        Real default_equilibrium_distance;
        Real native_contact_cutoff;
    };
}