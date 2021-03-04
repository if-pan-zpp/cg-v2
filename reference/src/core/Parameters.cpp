#include "core/Parameters.hpp"
#include "utils/Units.hpp"
using namespace CG;
using namespace std;

/* A helper for passing AminoAcid. */
static pair<AminoAcid, string> p(string const& amino_acid, string const& name) {
    return make_pair(AminoAcid(amino_acid), name);
}

Parameters::Parameters() {
    default_equilibrium_distance = 3.8 * angstrom;
    native_contact_cutoff = 0.0;

    for (auto const& amino_acid: AminoAcid::all_names) {
        atom_radii[p(amino_acid, "N")] = 1.64;
        atom_radii[p(amino_acid, "CA")] = 1.88;
        atom_radii[p(amino_acid, "C")] = 1.61;
        atom_radii[p(amino_acid, "O")] = 1.42;
    }

    /* TODO: Ask whether these are correct.\
     * Also, after the setting of the radii in cg.f there is some
     * sort of "correction" that I don't quite understand what it's for. */

    atom_radii[p("PRO", "CB")] = 1.88;
    atom_radii[p("PRO", "CG")] = 1.88;
    atom_radii[p("PRO", "CD")] = 1.88;

    atom_radii[p("GLN", "CB")] = 1.88;
    atom_radii[p("GLN", "CG")] = 1.88;
    atom_radii[p("GLN", "CD")] = 1.61;
    atom_radii[p("GLN", "OE1")] = 1.42;
    atom_radii[p("GLN", "NE2")] = 1.64;

    atom_radii[p("CYS", "CB")] = 1.88;
    atom_radii[p("CYS", "SG")] = 1.77;

    atom_radii[p("VAL", "CB")] = 1.88;
    atom_radii[p("VAL", "CG1")] = 1.88;
    atom_radii[p("VAL", "CG2")] = 1.88;

    atom_radii[p("PHE", "CB")] = 1.88;
    atom_radii[p("PHE", "CG")] = 1.88;
    atom_radii[p("PHE", "CD1")] = 1.61;
    atom_radii[p("PHE", "CD2")] = 1.76;
    atom_radii[p("PHE", "CE1")] = 1.76;
    atom_radii[p("PHE", "CE2")] = 1.76;
    atom_radii[p("PHE", "CZ")] = 1.76;

    atom_radii[p("MET", "CB")] = 1.88;
    atom_radii[p("MET", "CG")] = 1.88;
    atom_radii[p("MET", "SD")] = 1.77;
    atom_radii[p("MET", "CE")] = 1.88;

    atom_radii[p("ILE", "CB")] = 1.88;
    atom_radii[p("ILE", "CG1")] = 1.88;
    atom_radii[p("ILE", "CG2")] = 1.88;
    atom_radii[p("ILE", "CD1")] = 1.88;

    atom_radii[p("ASP", "CB")] = 1.88;
    atom_radii[p("ASP", "CG")] = 1.61;
    atom_radii[p("ASP", "OD1")] = 1.46;
    atom_radii[p("ASP", "OD1")] = 1.42;

    atom_radii[p("GLU", "CB")] = 1.88;
    atom_radii[p("GLU", "CG")] = 1.88;
    atom_radii[p("GLU", "CD")] = 1.61;
    atom_radii[p("GLU", "OE1")] = 1.46;
    atom_radii[p("GLU", "OE2")] = 1.42;

    atom_radii[p("LYS", "CB")] = 1.88;
    atom_radii[p("LYS", "CG")] = 1.88;
    atom_radii[p("LYS", "CD")] = 1.88;
    atom_radii[p("LYS", "CE")] = 1.88;
    atom_radii[p("LYS", "NZ")] = 1.64;

    atom_radii[p("ARG", "CB")] = 1.88;
    atom_radii[p("ARG", "CG")] = 1.88;
    atom_radii[p("ARG", "CD")] = 1.88;
    atom_radii[p("ARG", "NE")] = 1.64;
    atom_radii[p("ARG", "CZ")] = 1.61;
    atom_radii[p("ARG", "NH1")] = 1.64;
    atom_radii[p("ARG", "NH2")] = 1.64;

    atom_radii[p("SER", "CB")] = 1.88;
    atom_radii[p("SER", "OG")] = 1.46;

    atom_radii[p("THR", "CB")] = 1.88;
    atom_radii[p("THR", "OG1")] = 1.46;
    atom_radii[p("THR", "CG2")] = 1.88;

    atom_radii[p("TYR", "CB")] = 1.88;
    atom_radii[p("TYR", "CG")] = 1.61;
    atom_radii[p("TYR", "CD1")] = 1.76;
    atom_radii[p("TYR", "CD2")] = 1.76;
    atom_radii[p("TYR", "CE1")] = 1.76;
    atom_radii[p("TYR", "CE2")] = 1.76;
    atom_radii[p("TYR", "CZ")] = 1.61;
    atom_radii[p("TYR", "OH")] = 1.46;

    atom_radii[p("HIS", "CB")] = 1.88;
    atom_radii[p("HIS", "CG")] = 1.61;
    atom_radii[p("HIS", "ND1")] = 1.64;
    atom_radii[p("HIS", "CD2")] = 1.76;
    atom_radii[p("HIS", "CE1")] = 1.76;
    atom_radii[p("HIS", "NE2")] = 1.64;

    atom_radii[p("ASN", "CB")] = 1.88;
    atom_radii[p("ASN", "CG")] = 1.61;
    atom_radii[p("ASN", "OD1")] = 1.42;
    atom_radii[p("ASN", "ND2")] = 1.64;

    atom_radii[p("TRP", "CB")] = 1.88;
    atom_radii[p("TRP", "CG")] = 1.61;
    atom_radii[p("TRP", "CD1")] = 1.76;
    atom_radii[p("TRP", "NE1")] = 1.61;
    atom_radii[p("TRP", "CE2")] = 1.64;
    atom_radii[p("TRP", "CD2")] = 1.61;
    atom_radii[p("TRP", "CE3")] = 1.76;
    atom_radii[p("TRP", "CZ3")] = 1.76;
    atom_radii[p("TRP", "CZ2")] = 1.76;
    atom_radii[p("TRP", "CH2")] = 1.76;

    atom_radii[p("ALA", "CB")] = 1.88;

    atom_radii[p("LEU", "CB")] = 1.88;
    atom_radii[p("LEU", "CG")] = 1.88;
    atom_radii[p("LEU", "CD1")] = 1.88;
    atom_radii[p("LEU", "CD2")] = 1.88;

    /* Also, since the values are in multiples of 5 A, we need to
     * adjust these values. */
    for (auto& [name, radius]: atom_radii) {
        radius *= cg_f_unit;
    }
}
