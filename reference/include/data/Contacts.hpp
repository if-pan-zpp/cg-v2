#pragma once
#include "data/Primitives.hpp"
#include <vector>

namespace cg {
    class Contacts {
    public:
        /* We use pointers because rearranging a pair table when contacts
         * are added/removed is presumably expensive */
        std::vector<Pairs*> contacts;
    };
}