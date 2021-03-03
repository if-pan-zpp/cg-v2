#pragma once
#include "utils/Types.hpp"

namespace CG {
    /* This class provides an interface for the forces for use by the
     * integrators. */
    class Force {
    public:
        virtual ~Force() = 0;

        /* Compute forces and energies, or rather add them to the provided
         * arrays; if the pointer is null, force kernel should not compute
         * that quantity. */
        virtual void compute(Real3List *force, Real3List *energy) = 0;
    };
}