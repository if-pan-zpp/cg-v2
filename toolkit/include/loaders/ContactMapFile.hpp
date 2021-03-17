#pragma once
#include <istream>
#include <filesystem>
#include "model/NativeStructure.hpp"

namespace cg::toolkit {
    class ContactMapFile {
    public:
        ContactMapFile() = default;
        explicit ContactMapFile(std::istream& file);
        explicit ContactMapFile(std::filesystem::path const& path);
        NativeStructure ns;
    };
}