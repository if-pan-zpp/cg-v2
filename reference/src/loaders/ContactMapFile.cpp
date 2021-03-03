#include "loaders/ContactMapFile.hpp"
#include "utils/Units.hpp"
using namespace CG;
using namespace std;

ContactMapFile::ContactMapFile(istream &file) {
    int ncontacts, nresidues;
    file >> ns.offset >> nresidues >> ncontacts;

    ns.contacts = vector<NativeStructure::Contact>(ncontacts);
    ns.bond = RealList(nresidues);
    ns.dihedral = RealList(nresidues);

    for (Index i = 0; i < ncontacts; ++i) {
        /* Residue indices are given 1-indexed; we convert them to 0-indexed.
         * Also, we add the necessary offset. */
        auto& [res1, res2] = ns.contacts[i].residues;
        file >> res1;
        res1 -= 1;
        file >> res2;
        res2 -= 1;

        /* Bond distance is given in 5 Angstroms, per README.txt */
        auto& dist = ns.contacts[i].distance;
        file >> dist;
        dist *= 5.0 * angstrom;
    }

    for (Index i = 0; i < nresidues; ++i) {
        /* Angles are given in radians, per README.txt */
        file >> ns.bond[i];
        ns.bond[i] *= radian;
        file >> ns.dihedral[i];
        ns.dihedral[i] *= radian;
    }
}
