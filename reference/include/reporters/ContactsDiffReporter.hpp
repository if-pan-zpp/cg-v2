#pragma once
#include "reporters/Reporter.hpp"
#include "data/Results.hpp"
#include <map>

namespace cg::reference {
    class ContactsDiffReporter: public Reporter {
    public:
        ContactsDiffReporter(Results const &results, string filePath);

        void report(int step) override;
    private:
        Results const &results;

        map<unsigned, vector<QuasiAdiabatic::Contact>> refContacts;
    };
}
