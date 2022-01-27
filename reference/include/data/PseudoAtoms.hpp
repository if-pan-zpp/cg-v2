#pragma once
#include "data/Primitives.hpp"
#include "utils/RNG.hpp"
#include <unordered_map>

namespace cg::reference {   
    enum AACode : uint8_t;

    class PseudoAtoms {
    public:
        size_t n = 0; // number of pseudoatoms
        Reals3 pos;   // positions of pseudoatoms
        Reals3 vel;   // velocities of pseudoatoms multiplied by delta
        Reals mass, mass_inv, charge;
        Integers chainId;

        std::vector<AACode> aminoAcidCode;

        std::vector<std::string> type;
        std::unordered_map<std::string, std::pair<int, int>> typeRanges;

        // This method initializes velocities of residues, so that
        // the momentum is 0 and the average kinetic energy
        // corresponds to temperature (given in eps/kB units).
        void initMovement(toolkit::RNG &rng, Real temperature, Real delta);
    };

    enum AACode : uint8_t {
        ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY, HIS, ILE,
        LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL,
        NONE
    };
    const int NUM_AMINO_ACIDS = AACode::NONE;
    AACode aaCodeFromName(string const &name);
}
