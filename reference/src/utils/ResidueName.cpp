#include "utils/ResidueName.hpp"

using namespace CG;
using namespace std;

ResidueName::ResidueName(char c) {
    name_to_code = {
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
    code_to_name = {
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
    legacy_order = "GPQCASVTILNDKEMHFRYW";

    code = c;
}

ResidueName::ResidueName(const string &s) {
    code = name_to_code[s];
}

ResidueName::ResidueName(int i) {
    code = legacy_order[i];
}

ResidueName::operator char() const {
    return code;
}

ResidueName::operator std::string() const {
    return code_to_name[code];
}

ResidueName::operator int() const {
    return legacy_order.find(code);
}
