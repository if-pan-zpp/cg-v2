#pragma once
#include "utils/Types.hpp"
#include "model/Model.hpp"
#include "data/Parameters.hpp"
#include <unordered_map>
#include <unordered_set>
using namespace std;

namespace cg::toolkit {
    /* A "full" model. The necessity for this class comes from the fact that
     * deriving inter-residue contacts from full atomic data requires, well,
     * all atoms, and that we may wish to create a model from many PDB files,
     * derive contacts and then reduce FullModel to a normal Model. */
    class FullModel {
    public:
        FullModel() = default;

        /* Residues (with atoms) and chains. */
        struct Residue {
            unordered_map<string, Real3> atoms;
            AminoAcid type;
        };
        using Chain = vector<Residue>;
        unordered_map<Index, Chain> chains;

        /* Contacts (normal, SS bonds etc.). */
        using ResidueID = pair<Index, Index>;
        struct Contact {
            ResidueID res1, res2;
            Real distance;
            string type;
        };
        vector<Contact> contacts;
        using IndexPair = pair<Index, Index>;
        unordered_map<Index, unordered_set<IndexPair>> intraChainContacts;

        /* Operators for a disjoint union of full models. */
        FullModel &operator+=(FullModel const &fullModel2);
        FullModel operator+(FullModel const &fullModel2) const;

        /* Apply an affine transform. */
        void apply(RealAffine3 const &aff);

        /* Reduce into a Model. */
        Model reduce();

        /* Derive a contact map from full atomic data. */
        void deriveContactsFromAllAtoms(Parameters const &parameters);

    private:
        /* Reduce into a CG::Chain. */
        cg::toolkit::Chain reduceChain(Chain const &chain,
                                       unordered_set<IndexPair> const &intraContacts) const;
    };
}
