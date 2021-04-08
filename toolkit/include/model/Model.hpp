#pragma once
#include "model/Chain.hpp"
#include "data/Parameters.hpp"
#include <unordered_map>

namespace cg::toolkit {
    /* This class represents a molecular model (residues, SS bonds, angles etc.)
     * in an "intermediate" fashion, i.e. before being passed to the
     * System proper. */
    class Model {
    public:
        Model() = default;

        /* Construct a model from a single chain. */
        explicit Model(Chain const &chain);

        /* Model = a set of chains plus (possibly) extra structure.
         * We use a map in order to easily add/remove chains.
         * Note: not sure chaining with numbers is good or not. */
        std::unordered_map<Index, Chain> chains;

        /* Extra contacts (for example inter-chain bonds or SS bonds). */
        using ResidueID = std::pair<Index, Index>;
        struct Contact {
            ResidueID res1, res2;
            Real distance;
            std::string type;
        };
        std::vector<Contact> contacts;

        /* Operators for a disjoint sum of models. */
        Model &operator+=(Model const &model2);
        Model operator+(Model const &model2) const;

        /* Apply an affine transform to the entire model. */
        void apply(RealAffine3 const &aff);

        /* Derive a contact map from only CA atom positions. */
        void deriveContactsFromCaAtoms(Parameters const &parameters);
    };
}
