#pragma once
#include <istream>
#include <filesystem>
#include "data/Parameters.hpp"

namespace cg::toolkit {
    class ParameterFile {
    public:
        explicit ParameterFile(std::istream &file, bool loadMjMatrix = false);
        explicit ParameterFile(std::filesystem::path const &path, bool loadMjMatrix = false);
        Parameters parameters;
    };
}
