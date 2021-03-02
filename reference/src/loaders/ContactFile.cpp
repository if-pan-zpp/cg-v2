#include "loaders/ContactFile.hpp"
#include "math/Units.hpp"
using namespace CG;
using namespace std;

ContactFile::ContactFile(std::istream &file) {
    int ncontacts, nresidues;
    file >> offset >> nresidues >> ncontacts;

    contacts = vector<Contact>(ncontacts);
    angles = vector<Angles>(nresidues);

    for (int i = 0; i < ncontacts; ++i) {
        /* Residue indices are given 1-indexed; we convert them to 0-indexed. */
        auto& [res1, res2] = contacts[i].residues;
        file >> res1;
        --res1;
        file >> res2;
        --res2;

        /* Bond distance is given in 5 Angstrom multiples. */
        auto& dist = contacts[i].bond_distance;
        file >> dist;
        dist *= 5.0 * angstrom;
    }

    for (int i = 0; i < nresidues; ++i) {
        file >> angles[i].bond;
        angles[i].bond *= radian;
        file >> angles[i].dihedral;
        angles[i].dihedral *= radian;
    }
}
