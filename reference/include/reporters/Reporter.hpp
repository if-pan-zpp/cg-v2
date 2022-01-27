#pragma once

namespace cg::reference {
    class Reporter {
    public:
        virtual void report(int step) = 0;
    };
}
