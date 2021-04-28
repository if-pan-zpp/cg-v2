#include "data/Results.hpp"
using namespace cg::reference;

Results::Results() {
    clear();
}

void Results::clear() {
    potEnergy = 0.0;
    activeContacts = 0;
    qaContacts = NULL;
}
