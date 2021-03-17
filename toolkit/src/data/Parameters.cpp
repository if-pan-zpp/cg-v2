#include "data/Parameters.hpp"
#include "utils/Units.hpp"
using namespace cg::toolkit;
using namespace std;

/* A helper for passing AminoAcid. */
static pair<AminoAcid, string> p(string const& aminoAcid, string const& name) {
    return make_pair(AminoAcid(aminoAcid), name);
}

Parameters::Parameters() {
    defaultEquilibriumDistance = 3.8 * angstrom;
    nativeContactCutoff = 0.0;

    for (auto const& aminoAcid: AminoAcid::allNames) {
        atomRadii[p(aminoAcid, "N")] = 1.64;
        atomRadii[p(aminoAcid, "CA")] = 1.88;
        atomRadii[p(aminoAcid, "C")] = 1.61;
        atomRadii[p(aminoAcid, "O")] = 1.42;
    }

    /* TODO: Ask whether these are correct.\
     * Also, after the setting of the radii in cg.f there is some
     * sort of "correction" that I don't quite understand what it's for. */

    atomRadii[p("PRO", "CB")] = 1.88;
    atomRadii[p("PRO", "CG")] = 1.88;
    atomRadii[p("PRO", "CD")] = 1.88;

    atomRadii[p("GLN", "CB")] = 1.88;
    atomRadii[p("GLN", "CG")] = 1.88;
    atomRadii[p("GLN", "CD")] = 1.61;
    atomRadii[p("GLN", "OE1")] = 1.42;
    atomRadii[p("GLN", "NE2")] = 1.64;

    atomRadii[p("CYS", "CB")] = 1.88;
    atomRadii[p("CYS", "SG")] = 1.77;

    atomRadii[p("VAL", "CB")] = 1.88;
    atomRadii[p("VAL", "CG1")] = 1.88;
    atomRadii[p("VAL", "CG2")] = 1.88;

    atomRadii[p("PHE", "CB")] = 1.88;
    atomRadii[p("PHE", "CG")] = 1.88;
    atomRadii[p("PHE", "CD1")] = 1.61;
    atomRadii[p("PHE", "CD2")] = 1.76;
    atomRadii[p("PHE", "CE1")] = 1.76;
    atomRadii[p("PHE", "CE2")] = 1.76;
    atomRadii[p("PHE", "CZ")] = 1.76;

    atomRadii[p("MET", "CB")] = 1.88;
    atomRadii[p("MET", "CG")] = 1.88;
    atomRadii[p("MET", "SD")] = 1.77;
    atomRadii[p("MET", "CE")] = 1.88;

    atomRadii[p("ILE", "CB")] = 1.88;
    atomRadii[p("ILE", "CG1")] = 1.88;
    atomRadii[p("ILE", "CG2")] = 1.88;
    atomRadii[p("ILE", "CD1")] = 1.88;

    atomRadii[p("ASP", "CB")] = 1.88;
    atomRadii[p("ASP", "CG")] = 1.61;
    atomRadii[p("ASP", "OD1")] = 1.46;
    atomRadii[p("ASP", "OD2")] = 1.42;

    atomRadii[p("GLU", "CB")] = 1.88;
    atomRadii[p("GLU", "CG")] = 1.88;
    atomRadii[p("GLU", "CD")] = 1.61;
    atomRadii[p("GLU", "OE1")] = 1.46;
    atomRadii[p("GLU", "OE2")] = 1.42;

    atomRadii[p("LYS", "CB")] = 1.88;
    atomRadii[p("LYS", "CG")] = 1.88;
    atomRadii[p("LYS", "CD")] = 1.88;
    atomRadii[p("LYS", "CE")] = 1.88;
    atomRadii[p("LYS", "NZ")] = 1.64;

    atomRadii[p("ARG", "CB")] = 1.88;
    atomRadii[p("ARG", "CG")] = 1.88;
    atomRadii[p("ARG", "CD")] = 1.88;
    atomRadii[p("ARG", "NE")] = 1.64;
    atomRadii[p("ARG", "CZ")] = 1.61;
    atomRadii[p("ARG", "NH1")] = 1.64;
    atomRadii[p("ARG", "NH2")] = 1.64;

    atomRadii[p("SER", "CB")] = 1.88;
    atomRadii[p("SER", "OG")] = 1.46;

    atomRadii[p("THR", "CB")] = 1.88;
    atomRadii[p("THR", "OG1")] = 1.46;
    atomRadii[p("THR", "CG2")] = 1.88;

    atomRadii[p("TYR", "CB")] = 1.88;
    atomRadii[p("TYR", "CG")] = 1.61;
    atomRadii[p("TYR", "CD1")] = 1.76;
    atomRadii[p("TYR", "CD2")] = 1.76;
    atomRadii[p("TYR", "CE1")] = 1.76;
    atomRadii[p("TYR", "CE2")] = 1.76;
    atomRadii[p("TYR", "CZ")] = 1.61;
    atomRadii[p("TYR", "OH")] = 1.46;

    atomRadii[p("HIS", "CB")] = 1.88;
    atomRadii[p("HIS", "CG")] = 1.61;
    atomRadii[p("HIS", "ND1")] = 1.64;
    atomRadii[p("HIS", "CD2")] = 1.76;
    atomRadii[p("HIS", "CE1")] = 1.76;
    atomRadii[p("HIS", "NE2")] = 1.64;

    atomRadii[p("ASN", "CB")] = 1.88;
    atomRadii[p("ASN", "CG")] = 1.61;
    atomRadii[p("ASN", "OD1")] = 1.42;
    atomRadii[p("ASN", "ND2")] = 1.64;

    atomRadii[p("TRP", "CB")] = 1.88;
    atomRadii[p("TRP", "CG")] = 1.61;
    atomRadii[p("TRP", "CD1")] = 1.76;
    atomRadii[p("TRP", "NE1")] = 1.61;
    atomRadii[p("TRP", "CE2")] = 1.64;
    atomRadii[p("TRP", "CD2")] = 1.61;
    atomRadii[p("TRP", "CE3")] = 1.76;
    atomRadii[p("TRP", "CZ3")] = 1.76;
    atomRadii[p("TRP", "CZ2")] = 1.76;
    atomRadii[p("TRP", "CH2")] = 1.76;

    atomRadii[p("ALA", "CB")] = 1.88;

    atomRadii[p("LEU", "CB")] = 1.88;
    atomRadii[p("LEU", "CG")] = 1.88;
    atomRadii[p("LEU", "CD1")] = 1.88;
    atomRadii[p("LEU", "CD2")] = 1.88;

    /* Also, since the values are in multiples of 5 A, we need to
     * adjust these values. */
    for (auto& [name, radius]: atomRadii) {
        radius *= f77unit;
    }
}
