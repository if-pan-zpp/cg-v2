#pragma once
#include <iostream>
#include <mdk/tools/param/Parameters.hpp>

namespace mdk::param {
    class LegacyParser {
    public:
        Parameters read(std::istream& is);
        std::ostream& write(std::ostream& os, Parameters const& data);
    };
}