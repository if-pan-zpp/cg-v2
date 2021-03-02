#include "loaders/PDBFile.hpp"
#include "math/Types.hpp"
#include "math/Units.hpp"
#include <stdexcept>
#include <map>

using namespace std;
using namespace loaders;
using namespace math;

/* Yield a view of characters i..j (1-indexed) of s.
 * Note: 1-indexing is due to the specs in PDB file format
 * reference. */
static string view(string const &s, unsigned i, unsigned j) {
    if (i > s.size() || j > s.size() || i > j)
        throw runtime_error("invalid view requested");

    return s.substr(i - 1, j - i + 1);
}

/* Much the same, but a single char. Mostly for aesthetics. */
static char view(string const &s, unsigned i) {
    if (i > s.size())
        throw runtime_error("invalid view requested");

    return s[i - 1];
}

/* Remove whitespace from the left. */
static inline string ltrim(string s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), [](unsigned char c) -> auto {
        return !std::isspace(c);
    }));
    return s;
}

/* Remove whitespace from the right. */
static string rtrim(string s) {
    s.erase(find_if(s.rbegin(), s.rend(), [](unsigned char c) -> auto {
        return !std::isspace(c);
    }).base(), s.end());
    return s;
}

/* Remove whitespace from both ends. */
static string trim(string s) {
    return rtrim(ltrim(move(s)));
}

/* Map from residue names to codes. */
static map<string, char> residue_name_to_code = {
    {"ALA", 'A'},
    {"ARG", 'R'},
    {"ASN", 'N'},
    {"ASP", 'D'},
    {"CYS", 'C'},
    {"GLU", 'E'},
    {"GLN", 'Q'},
    {"GLY", 'G'},
    {"HIS", 'H'},
    {"ILE", 'I'},
    {"LEU", 'L'},
    {"LYS", 'K'},
    {"MET", 'M'},
    {"PHE", 'F'},
    {"PRO", 'P'},
    {"SER", 'S'},
    {"THR", 'T'},
    {"TRP", 'W'},
    {"TYR", 'Y'},
    {"VAL", 'V'}
};

PDBFile::PDBFile(istream &file, bool unwrap) {
    /* Note: we don't parse TER records. */
    /* TODO: lunwrap */

    string line;
    while (getline(file, line)) {
        /* Instead of requiring exactly 80 chars, we just ignore the rest. */
        if (line.size() < 80)
            throw runtime_error("PDB - line size < 80");

        auto record_name = trim(view(line, 1, 6));
        if (record_name == "ATOM") {
            /* Parse the necessary columns.
             * Note: Residue seq nums are 1-indexed by default, we change it. */
            string atom_name = trim(view(line, 13, 16));
            char alt_location = view(line, 17);
            string residue_name = trim(view(line, 18, 20));
            char chain_id = view(line, 22);
            int residue_seq_num = stoi(view(line, 23, 26)) - 1;
            Real x = stod(view(line, 31, 38)) * angstrom;
            Real y = stod(view(line, 39, 46)) * angstrom;
            Real z = stod(view(line, 47, 54)) * angstrom;

            /* cg.f:5080 - Essentially we ignore alternative locations, if any
             * are present. */
            if (alt_location != 'A' && alt_location != ' ')
                continue;

            /* TODO: cg.f:5061-5064 describes "filtering out DNA chains" but I
             * couldn't find anything resembling DNA chains in provided PDB
             * files, and the code is strange. */

            /* Locate the chain; if such chain doesn't exist, operator[]
             * creates it. */
            auto& chain = chains[chain_id];

            /* Locate the residue. */
            if (residue_seq_num > chain.residues.size())
                throw runtime_error("PDB - incorrect residue order");
            if (residue_seq_num == chain.residues.size())
                chain.residues.emplace_back();
            auto& residue = chain.residues[residue_seq_num];

            /* Insert the atom. */
            auto code_iter = residue_name_to_code.find(residue_name);
            if (code_iter == residue_name_to_code.end())
                throw runtime_error("PDB - incorrect residue name");
            residue.residue_code = code_iter->second;
            residue.atoms[atom_name] = Real3(x, y, z);
        }
        else if (record_name == "SSBOND") {
            /* Parse relevant fields
             * Note: Residue seq nums are 1-indexed by default, we change it. */
            char cys1_chain_id = view(line, 16);
            int cys1_residue_seq_num = stoi(view(line, 18, 21)) - 1;
            char cys2_chain_id = view(line, 30);
            int cys2_residue_seq_num = stoi(view(line, 32, 35)) - 1;
            Real bond_distance = stod(view(line, 74, 78)) * angstrom;

            /* Insert SS bond */
            ssbonds.emplace_back((SSBond) {
                .cys1 = {cys1_chain_id, cys1_residue_seq_num},
                .cys2 = {cys2_chain_id, cys2_residue_seq_num},
                .bond_distance = bond_distance
            });
        }
        else if (record_name == "CRYST1" && !cryst1) {
            /* Parse relevant fields. cg.f:5066-5076, we only want the sides
             * and ignore the angles and space group and Z value. */
            Real a = stod(view(line, 7, 15)) * angstrom;
            Real b = stod(view(line, 16, 24)) * angstrom;
            Real c = stod(view(line, 25, 33)) * angstrom;

            cryst1 = (Cryst1) {.size = Real3(a, b, c)};
        }
        else if (record_name == "END") {
            break;
        }
    }
}