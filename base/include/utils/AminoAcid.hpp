#pragma once
#include <vector>
#include <string>
#include <unordered_map>

namespace cg {
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

        static std::string allCodes;
        static std::vector<std::string> allNames;
        static int numAminoAcids;

        bool operator==(AminoAcid const& aminoAcid2) const;

    private:
        static std::unordered_map<std::string, char> nameToCode;
        static std::unordered_map<char, std::string> codeToName;

        char code = 0;
        std::string name = "";
    };
}

namespace std {
    template<>
    struct hash<cg::AminoAcid> {
        size_t operator()(cg::AminoAcid const& aminoAcid) const {
            return std::hash<char>()((char)aminoAcid);
        }
    };
}