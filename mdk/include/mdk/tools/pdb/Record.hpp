#pragma once
#include <Eigen/Core>
#include <string>
#include <variant>

namespace mdk::pdb::records {
    class Atom {
    public:
        int serialNum = 0;
        std::string atomName;
        char altLocation = ' ';
        std::string residueName;
        char chainID = 'A';
        int residueSeqNum = 0;
        char insertionCode = ' ';
        Eigen::Vector3d pos;
        double occupancy = 1.0;
        double tempFactor = 0.0;
        std::string element;
        std::string charge;
    };

    class SSBond {
    public:
        int serialNum = 0;

        struct PerResidueData {
            char chainId = 'A';
            int residueSeqNum = 0;
            char insertionCode = ' ';
            std::string symmetryOp;
        };
        PerResidueData res[2];

        double dist0 = 0.0;
    };

    class Cryst1 {
    public:
        Eigen::Vector3d size;
        Eigen::Vector3d angles;
        std::string spaceGroup;
        int z = 0;
    };

    class Ter {
    public:
        int serialNum = 0;
        std::string residueName;
        char chainId = 'A';
        int residueSeqNum = 0;
        char insertionCode = ' ';
    };

    class Model {
    public:
        int serialNum = 0;
    };

    class End {
    public:
    };

    using Record = std::variant<std::monostate,
        Atom, SSBond, Cryst1, End, Model, Ter>;
}