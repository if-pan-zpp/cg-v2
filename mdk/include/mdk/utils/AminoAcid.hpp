#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace mdk {
    class AminoAcid {
    public:
        AminoAcid();
        AminoAcid(char c);
        AminoAcid(std::string const &s);

        static bool isAminoAcid(std::string const& s);

        explicit operator char() const;
        explicit operator std::string() const;

        static std::unordered_map<std::string, char> nameToCode;
        static std::unordered_map<char, std::string> codeToName;
        static std::unordered_set<AminoAcid> allAminoAcids;
        static std::unordered_set<std::string> allNames;
        static constexpr int numAminoAcids = 20;

        bool operator==(AminoAcid const &aminoAcid2) const;

    private:
        char code;
        std::string name;
    };
}

namespace std {
    template<>
    struct hash<mdk::AminoAcid> {
        size_t operator()(mdk::AminoAcid const &aminoAcid) const {
            return std::hash<char>()((char)aminoAcid);
        }
    };
}
