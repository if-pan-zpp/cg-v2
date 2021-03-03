#pragma once
#include "utils/Types.hpp"
#include "model/Model.hpp"
#include "model/NativeStructure.hpp"
#include <istream>
#include <vector>
#include <string>
#include <unordered_map>
#include <optional>

namespace CG {
    /* PDB file loader */
    class PDBFile {
    public:
        /* Read a PDB file from a given stream. If relevant CRYST1
         * field is present, the crystallographic structure may be
         * unwrapped. */
        explicit PDBFile(std::istream& file, bool unwrap = false);

        /* Construct a model from the PDB file. */
        Model model(bool from_all_atoms = false);

    private:
        struct Residue {
            /* We store for each atom only its position */
            std::unordered_map<std::string, Real3> atoms;
            char residue_code;
        };
        struct Chain {
            std::vector<Residue> residues;
        };
        std::unordered_map<char, Chain> chains;

        struct SSBond {
            /* Chain identifier and residue index. */
            std::pair<char, Index> cys1, cys2;
            Real bond_distance;
        };
        std::vector<SSBond> ssbonds;

        /* Crystallographic data */
        struct Cryst1 {
            Real3 size;
        };
        std::optional<Cryst1> cryst1;
    };
}