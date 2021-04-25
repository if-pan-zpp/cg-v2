#pragma once
#include <vector>
#include <string>
#include <Eigen/Core>
#include <mdk/tools/model/CoarseModel.hpp>
#include <optional>

namespace mdk {
    class AtomicModel {
    public:
        struct Atom;
        struct Residue;
        struct Chain;

        struct Atom {
            int serial, idxInRes;
            Residue *res;
            Eigen::Vector3d pos;
            std::string type;
        };
        std::unordered_map<int, Atom> atoms;

        Atom& addAtom(int idx, Residue *res = nullptr);

        struct Residue {
            int serial, idxInChain;
            Chain *chain;
            std::string type;
            std::vector<Atom*> atoms;

            Atom* find(std::string const& type);
        };
        std::unordered_map<int, Residue> residues;

        Residue& addResidue(int idx, Chain *chain = nullptr);

        struct Chain {
            char serial;
            std::vector<Residue*> residues;
        };
        std::unordered_map<char, Chain> chains;

        Chain& addChain(char idx);

        struct Contact {
            int idx;
            Atom* atom[2];
            std::string type;
            double dist0;
        };
        std::vector<Contact> contacts;

        Contact& addContact();

        std::optional<Eigen::Vector3d> cell;

    public:
        AtomicModel() = default;
        explicit AtomicModel(CoarseModel const& coarse);

        void addNativeContacts(bool onlyCA = false);
        CoarseModel coarsen();
    };
}
