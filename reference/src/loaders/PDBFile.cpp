#include "loaders/PDBFile.hpp"
#include "utils/Types.hpp"
#include "utils/Units.hpp"
#include "utils/AminoAcid.hpp"
#include <stdexcept>
#include <map>

using namespace std;
using namespace CG;

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

PDBFile::PDBFile(istream &file) {
    /* Note: we don't parse TER records. As for the magic numbers, see PDB
     * file reference.
     * TODO: implement unwrapping */

    /* Map from the chain id's within the file to the ones in the full_model. */
    std::unordered_map<char, Index> chain_id_map;

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
            if (chain_id_map.find(chain_id) == chain_id_map.end())
                chain_id_map[chain_id] = chain_id_map.size();
            auto &chain = full_model.chains[chain_id_map[chain_id]];

            /* Locate the residue. */
            if (residue_seq_num > chain.size())
                throw runtime_error("PDB - incorrect residue order");
            auto &residue = chain.emplace_back();

            /* Insert the atom. */
            residue.type = residue_name;
            residue.atoms[atom_name] = Real3(x, y, z);
        } else if (record_name == "SSBOND") {
            /* Parse relevant fields
             * Note: Residue seq nums are 1-indexed by default, we change it. */
            char cys1_chain_id = view(line, 16);
            int cys1_seq_num = stoi(view(line, 18, 21)) - 1;
            char cys2_chain_id = view(line, 30);
            int cys2_seq_num = stoi(view(line, 32, 35)) - 1;
            Real distance = stod(view(line, 74, 78)) * angstrom;

            /* Insert SS bond */
            full_model.contacts.push_back((FullModel::Contact) {
                .res1 = { chain_id_map[cys1_chain_id], cys1_seq_num },
                .res2 = { chain_id_map[cys2_chain_id], cys2_seq_num },
                .distance = distance,
                .type = "ssbond"
            });
        } else if (record_name == "CRYST1" && !cryst1) {
            /* Parse relevant fields. cg.f:5066-5076, we only want the sides
             * and ignore the angles and space group and Z value. */
            Real a = stod(view(line, 7, 15)) * angstrom;
            Real b = stod(view(line, 16, 24)) * angstrom;
            Real c = stod(view(line, 25, 33)) * angstrom;

            cryst1 = Real3(a, b, c);
        } else if (record_name == "END") {
            break;
        }
    }
}