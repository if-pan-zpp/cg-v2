#pragma once
#include "forces/Force.hpp"
#include "data/NativeStructure.hpp"
#include "data/Topology.hpp"

namespace cg::reference {
    /* CPC14.pdf, 4.1
     * Note: here we may forego neighborhood lists */
    class NativeContacts: public Force {
    private:
        NativeStructure const* ns;
        Topology const* top;

    public:
        Real depth = 1.0;

        NativeContacts(NativeStructure const &ns, Topology const &top);

        void compute(Real &energy, Reals3 &forces) override;
    };
}

