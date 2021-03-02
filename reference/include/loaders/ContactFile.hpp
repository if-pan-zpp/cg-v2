#pragma once
#include <istream>
#include <vector>
#include "math/Types.hpp"

namespace loaders {
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
            math::Real bond_distance;
        };
        std::vector<Contact> contacts;

        /* Bond, dihedral angles */
        struct Angles {
            math::Real bond;
            math::Real dihedral;
        };
        std::vector<Angles> angles;
    };
}