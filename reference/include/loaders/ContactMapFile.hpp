#pragma once
#include <istream>
#include "model/NativeStructure.hpp"

namespace CG {
    class ContactMapFile {
    public:
        explicit ContactMapFile(std::istream& file);
        NativeStructure ns;
    };
}