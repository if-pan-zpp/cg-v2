#pragma once
#include <istream>

namespace geometry {
    class PDBFile {
    public:
        PDBFile(std::istream& file);

    };
}