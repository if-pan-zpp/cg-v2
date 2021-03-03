#pragma once
#include <istream>
#include <vector>
#include "math/Vector.hpp"

namespace CG {
    class ContactFile {
    public:
        ContactFile() = default;
        explicit ContactFile(std::istream& file);

    private:
        /* "Reading frame offset" */
        int offset;

        /* Contacts in the map */
        struct Contact {
            std::pair<int, int> residues;
            Real bond_distance;
        };
        std::vector<Contact> contacts;

        /* Bond, dihedral angles */
        struct Angles {
            Real bond;
            Real dihedral;
        };
        std::vector<Angles> angles;
    };
}