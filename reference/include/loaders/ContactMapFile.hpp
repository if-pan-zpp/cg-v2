#pragma once
#include <istream>
#include "model/NativeStructure.hpp"

namespace CG {
    class ContactMapFile {
    public:
        ContactMapFile() = default;
        explicit ContactMapFile(std::istream& file);
        NativeStructure ns;
    };
}