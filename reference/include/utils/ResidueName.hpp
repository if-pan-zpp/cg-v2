#pragma once
#include <string>
#include <map>

namespace CG {
#define NUM_AMINOACIDS 20

    /* A generalized residue name class. */
    class ResidueName {
    public:
        explicit ResidueName(char c);
        explicit ResidueName(std::string const& s);
        explicit ResidueName(int i);

        explicit operator char() const;
        explicit operator std::string() const;
        explicit operator int() const;

    private:
        char code;

        static std::map<std::string, char> name_to_code;
        static std::map<char, std::string> code_to_name;
        static std::string legacy_order;
    };
}