#pragma once

#include "core/System.hpp"
#include <vector>
#include <map>

namespace CG {
    class Geometry;

    /* This class contains the particle pairs within a specified cutoff
     * distance (possibly with some spurious ones). */
    class Neighborhoods {
    public:
        IndexPairList pairs;

        /* Restrict the neighborhoods to a lesser cutoff and pad. */
//        Neighborhoods restrict(Real sub_cutoff, Real sub_pad) const;

        Neighborhoods() = default;

    private:
        friend class Geometry;

        /* Construct neighborhoods with a given cutoff and a "pad", which
         * is meant to buffer particle displacements and decrease the number
         * of times the list must be recomputed at the cost of having
         * more spurious pairs.
         *
         * More specifically, initially the list contains all the pairs within
         * the distance of cutoff+2*pad, and until the maximal relative
         * displacement of particles exceeds pad, all the pairs within the
         * distance of cutoff are guaranteed to be correct.
         *
         * include4 is lii4, essentially. */
        Neighborhoods(Geometry *geometry, Real cutoff, Real pad,
                bool include4 = true);

        Geometry *geometry;
        System *system;

        Real cutoff, pad;
        bool include4;
        Real max_correct_dist;

        void update_max_correct_dist();
        void update();

        /* Particle positions at the time of computing the list. */
        Real3List reference_pos;
        /* Box shape at the time of computing the list. */
        Real3 reference_box_shape;
    };

    /* This class is responsible for implementing the metric on the space
     * and for the construction of neighborhood lists.
     * Note: <Geometry> class does not contain the positions and velocities. */
    class Geometry {
    public:
        /* A geometry on a system. */
        explicit Geometry(System *system);

        /* Construct neighborhoods. The value being a pointer will allow us to
         * update the lists without having to reassign them. */
        const Neighborhoods *
        make_neighborhoods(Real cutoff, Real pad, bool include4 = true);

        /* Metric of the geometry. */
        [[nodiscard]] inline Real3 diff(Real3 p, Real3 q) const;

        /* Set/update the "simulation box" shape (i.e. the cell which is
         * repeated in space). The value of \leq 0 is interpreted as
         * infinity. */
        void set_simulation_box(Real3 box);

        /* Update the geometry (specifically the neighborhoods, if the particles
         * have been sufficiently displaced or due to the simulation box shape
         * changes. */
        void update();

    private:
        friend class Neighborhoods;

        System *system;

        /* All stored neighborhoods. Pointers to elements of a map shall not be
         * invalidated (unlike, for example, with std::vector). Also, we will
         * be able not to create extraneous neighborhoods. */
        using NeighborhoodsParams = std::tuple<Real, Real, bool>;
        std::map<NeighborhoodsParams, Neighborhoods> all_neighborhoods;

        /* PBC data. pbc = true iff the box is finite along the specified
         * axis. */
        std::array<bool, 3> pbc;
        Real3 box_shape, box_shape_inv;
    };
}