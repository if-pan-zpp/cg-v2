#include "forces/local/HarmonicTethers.hpp"
using namespace cg::reference;

HarmonicTethers::HarmonicTethers(PseudoAtoms const& _pseudoAtoms,
                                 NativeStructure const& _ns):
    pseudoAtoms(_pseudoAtoms),
    ns(_ns) {
    
}

void HarmonicTethers::compute(Real &energy, Reals3 &forces) {

}
