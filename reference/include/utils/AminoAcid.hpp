#pragma once
#include <vector>
#include <string>
#include <unordered_map>

namespace CG {
    /* Amino acid name -- used for representing full-name repr and code repr
     * uniformly. */
    class AminoAcid {
    public:
        /* We shall use implicit operators for brevity. */
        AminoAcid() = default;
        AminoAcid(char c);
        AminoAcid(std::string const& s);

        explicit operator char() const;
        explicit operator std::string() const;

        static std::string all_codes;
        static std::vector<std::string> all_names;
        static int num_amino_acids;

        bool operator==(AminoAcid const& amino_acid2) const;

    private:
        static std::unordered_map<std::string, char> name_to_code;
        static std::unordered_map<char, std::string> code_to_name;

        char code;
        std::string name;
    };
}

namespace std {
    template<>
    struct hash<CG::AminoAcid> {
        size_t operator()(CG::AminoAcid const& amino_acid) const {
            return std::hash<char>()((char)amino_acid);
        }
    };
}