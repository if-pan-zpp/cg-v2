#include "utils/AminoAcid.hpp"

using namespace cg::toolkit;
using namespace std;

string AminoAcid::allCodes = "GPQCASVTILNDKEMHFRYW";

vector<string> AminoAcid::allNames = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
};

unordered_map<char, string> AminoAcid::codeToName = {
        {'A', "ALA"},
        {'R', "ARG"},
        {'N', "ASN"},
        {'D', "ASP"},
        {'C', "CYS"},
        {'E', "GLU"},
        {'Q', "GLN"},
        {'G', "GLY"},
        {'H', "HIS"},
        {'I', "ILE"},
        {'L', "LEU"},
        {'K', "LYS"},
        {'M', "MET"},
        {'F', "PHE"},
        {'P', "PRO"},
        {'S', "SER"},
        {'T', "THR"},
        {'W', "TRP"},
        {'Y', "TYR"},
        {'V', "VAL"}
};

unordered_map<string, char> AminoAcid::nameToCode = {
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

int AminoAcid::numAminoAcids = 20;

AminoAcid::AminoAcid(char c) {
    code = c;
    name = codeToName[c];
}

AminoAcid::AminoAcid(const string &s) {
    name = s;
    code = nameToCode[s];
}

AminoAcid::operator char() const {
    return code;
}

AminoAcid::operator std::string() const {
    return name;
}

bool AminoAcid::operator==(AminoAcid const& aminoAcid2) const {
    return code == aminoAcid2.code;
}
