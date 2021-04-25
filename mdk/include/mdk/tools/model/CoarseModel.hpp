#pragma once
#include <vector>
#include <string>
#include <optional>
#include <Eigen/Core>
#include <mdk/utils/Units.hpp>
#include <mdk/utils/Random.hpp>
#include <mdk/tools/model/Topology.hpp>

namespace mdk {
    class CoarseModel {
    public:
        struct Residue;
        struct Chain;

        struct Residue {
            int idx, chainIdx;
            Eigen::Vector3d pos;
            std::string type;
        };
        std::vector<Residue> residues;
        Residue& addResidue(Chain *chain = nullptr);

        struct Chain {
            int idx;
            int first, last;
            std::vector<double> tethers;
            std::vector<int> structuredParts;
        };
        std::vector<Chain> chains;
        Chain& addChain();

        struct Contact {
            int idx;
            int res[2];
            std::string type;
            double dist0;
        };
        std::vector<Contact> contacts;
        Contact& addContact();

        struct StructuredPart {
            int idx;
            int off, len;
            std::vector<double> angle, dihedral;
        };
        std::vector<StructuredPart> structuredParts;
        StructuredPart& addSP();

    public:
        void morphIntoLine();
        std::optional<Topology> morphIntoSAW(Random& rand, bool withPBC, double density, double intersectionDist);

    private:
        std::vector<std::pair<Residue*, Residue*>> nonlocalPairs();
    };
}
