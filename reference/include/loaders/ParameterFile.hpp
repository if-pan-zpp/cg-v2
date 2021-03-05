#pragma once
#include <istream>
#include <filesystem>
#include "core/Parameters.hpp"

namespace CG {
    class ParameterFile {
    public:
        explicit ParameterFile(std::istream& file, bool load_mj_matrix = false);
        explicit ParameterFile(std::filesystem::path const& path, bool load_mj_matrix = false);
        Parameters parameters;
    };
}