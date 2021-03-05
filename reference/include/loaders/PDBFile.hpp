#pragma once
#include "utils/Types.hpp"
#include "model/FullModel.hpp"
#include <istream>
#include <vector>
#include <string>
#include <unordered_map>
#include <optional>
#include <filesystem>

namespace CG {
    /* PDB file loader */
    class PDBFile {
    public:
        /* Read a PDB file from a given stream. If relevant CRYST1
         * field is present.
         * Note: figure out what lunwrap does in cg.f */
        explicit PDBFile(std::istream& file);
        explicit PDBFile(std::filesystem::path const& path);

        FullModel full_model;
        std::optional<Real3> cryst1;
    };
}