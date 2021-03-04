#include "utils/AminoAcid.hpp"

using namespace CG;
using namespace std;

string AminoAcid::all_codes = "GPQCASVTILNDKEMHFRYW";

vector<string> AminoAcid::all_names = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
};

unordered_map<char, string> AminoAcid::code_to_name = {
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

unordered_map<string, char> AminoAcid::name_to_code = {
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

int AminoAcid::num_amino_acids = 20;

AminoAcid::AminoAcid(char c) {
    code = c;
    name = code_to_name[c];
}

AminoAcid::AminoAcid(const string &s) {
    name = s;
    code = name_to_code[s];
}

AminoAcid::operator char() const {
    return code;
}

AminoAcid::operator std::string() const {
    return name;
}

bool AminoAcid::operator==(AminoAcid const& amino_acid2) const {
    return code == amino_acid2.code;
}
