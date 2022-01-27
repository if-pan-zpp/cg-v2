#include "loaders/PDBFile.hpp"
#include "utils/Units.hpp"
#include <stdexcept>
#include <fstream>
#include <map>

using namespace cg::toolkit;
using namespace std;

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
    std::unordered_map<char, Index> chainIdMap;

    string line;
    while (getline(file, line)) {
        auto recordName = trim(view(line, 1, 6));
        if (recordName == "ATOM") {
            /* Parse the necessary columns.
             * Note: Residue seq nums are 1-indexed by default, we change it. */
            string atomName = trim(view(line, 13, 16));
            char altLocation = view(line, 17);
            string residueName = trim(view(line, 18, 20));
            char chainId = view(line, 22);
            int residueSeqNum = stoi(view(line, 23, 26)) - 1;
            Real x = stod(view(line, 31, 38)) * angstrom;
            Real y = stod(view(line, 39, 46)) * angstrom;
            Real z = stod(view(line, 47, 54)) * angstrom;

            /* cg.f:5080 - Essentially we ignore alternative locations, if any
             * are present. */
            if (altLocation != 'A' && altLocation != ' ')
                continue;

            /* TODO: cg.f:5061-5064 describes "filtering out DNA chains" but I
             * couldn't find anything resembling DNA chains in provided PDB
             * files, and the code is strange. */

            /* Locate the chain; if such chain doesn't exist, operator[]
             * creates it. */
            if (chainIdMap.find(chainId) == chainIdMap.end())
                chainIdMap[chainId] = chainIdMap.size();
            auto &chain = fullModel.chains[chainIdMap[chainId]];

            /* Locate the residue. */
            if (residueSeqNum > chain.size())
                throw runtime_error("PDB - incorrect residue order");
            if (residueSeqNum == chain.size())
                chain.emplace_back();
            auto &residue = chain[residueSeqNum];

            /* Insert the atom. */
            residue.type = residueName;
            residue.atoms[atomName] = Real3(x, y, z);
        } else if (recordName == "SSBOND") {
            /* Parse relevant fields
             * Note: Residue seq nums are 1-indexed by default, we change it. */
            char cys1ChainId = view(line, 16);
            int cys1SeqNum = stoi(view(line, 18, 21)) - 1;
            char cys2ChainId = view(line, 30);
            int cys2SeqNum = stoi(view(line, 32, 35)) - 1;
            Real distance = stod(view(line, 74, 78)) * angstrom;

            /* Insert SS bond */
            fullModel.contacts.push_back((FullModel::Contact) {
                .res1 = {chainIdMap[cys1ChainId], cys1SeqNum },
                .res2 = {chainIdMap[cys2ChainId], cys2SeqNum },
                .distance = distance,
                .type = "ssbond"
            });
        } else if (recordName == "CRYST1" && !cryst1) {
            /* Parse relevant fields. cg.f:5066-5076, we only want the sides
             * and ignore the angles and space group and Z value. */
            Real a = stod(view(line, 7, 15)) * angstrom;
            Real b = stod(view(line, 16, 24)) * angstrom;
            Real c = stod(view(line, 25, 33)) * angstrom;

            cryst1 = Real3(a, b, c);
        } else if (recordName == "END") {
            break;
        }
    }
}

PDBFile::PDBFile(const filesystem::path &path) {
    auto filestream = ifstream(path);
    *this = PDBFile(filestream);
}
