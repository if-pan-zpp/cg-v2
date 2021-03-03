#pragma once
#include <istream>
#include <unordered_map>
#include <string>
#include <filesystem>
#include "loaders/ContactMapFile.hpp"

namespace CG {
    /* Sequence file loader */
    class SequenceFile {
    public:
        /* This loader needs an explicit path, because we interpret
         * the contact map paths to be relative to the sequence file
         * and not the caller. */
        explicit SequenceFile(std::filesystem::path const& path);

    private:
        /* An array with contact files, so as to not load them
         * multiple times. */
        std::unordered_map<std::string, ContactMapFile> map_files;

        /* Chains. */
        struct Chain {
            std::string residue_codes;
            std::vector<ContactMapFile*> map_files;
        };
        std::vector<Chain> chains;
    };
}