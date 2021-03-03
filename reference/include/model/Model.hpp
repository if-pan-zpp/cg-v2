#pragma once
#include "model/Chain.hpp"
#include <unordered_map>

namespace CG {
    /* This class represents a molecular model (residues, SS bonds, angles etc.)
     * in an "intermediate" fashion, i.e. before being passed to the System
     * proper. */
    class Model {
    public:
        Model() = default;

        /* Construct a model from a single chain. */
        explicit Model(Chain const& chain, std::string const& name = "");

        /* Model = a set of chains plus (possibly) extra structure.
         * We use a map in order to easily add/remove chains. */
        std::unordered_map<Index, Chain> chains;

        struct Contact {
            std::pair<Index, Index> res1, res2;
            Real distance;
            std::string type;
        };
        std::vector<Contact> extra_contacts;

        /* Operators for a disjoint sum of models. */
        Model& operator+=(Model const& model2);
        Model operator+(Model const& model2) const;

        /* Apply an affine transform to the entire model. */
        void apply(RealAffine3 const& aff);
    };
}