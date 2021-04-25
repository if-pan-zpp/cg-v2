#pragma once
#include <vector>

namespace mdk::cmap {
    class ContactMap {
    public:
        int len, offset;
        struct Contact {
            int res[2];
            double dist0;
        };

        std::vector<Contact> contacts;
        std::vector<double> angle, dihedral;
    };
}