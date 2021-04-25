#include "utils/ResidueTypes.hpp"
#include "utils/Units.hpp"
using namespace mdk;

ResidueTypes mdk::createTypes() {
    ResidueTypes types;

    for (auto const& acid: AminoAcid::allAminoAcids) {
        auto& acidInfo = types[acid];
        acidInfo.atomInfo["N"].radius = 1.64;
        acidInfo.atomInfo["CA"].radius = 1.88;
        acidInfo.atomInfo["C"].radius = 1.61;
        acidInfo.atomInfo["O"].radius = 1.42;
        acidInfo.atomInfo["OXT"].radius = 1.42; // ??

        for (auto& [type, entry]: acidInfo.atomInfo) {
            entry.inBackbone = true;
        }
    }

    auto& PRO = types[(AminoAcid)"PRO"];
    PRO.atomInfo["CB"].radius = 1.88;
    PRO.atomInfo["CG"].radius = 1.88;
    PRO.atomInfo["CD"].radius = 1.88;

    auto& GLN = types[(AminoAcid)"GLN"];
    GLN.atomInfo["CB"].radius = 1.88;
    GLN.atomInfo["CB"].radius = 1.88;
    GLN.atomInfo["CG"].radius = 1.88;
    GLN.atomInfo["CD"].radius = 1.61;
    GLN.atomInfo["OE1"].radius = 1.42;
    GLN.atomInfo["NE2"].radius = 1.64;

    auto& CYS = types[(AminoAcid)"CYS"];
    CYS.atomInfo["CB"].radius = 1.88;
    CYS.atomInfo["SG"].radius = 1.77;

    auto& VAL = types[(AminoAcid)"VAL"];
    VAL.atomInfo["CB"].radius = 1.88;
    VAL.atomInfo["CG1"].radius = 1.88;
    VAL.atomInfo["CG2"].radius = 1.88;

    auto& PHE = types[(AminoAcid)"PHE"];
    PHE.atomInfo["CB"].radius = 1.88;
    PHE.atomInfo["CG"].radius = 1.88;
    PHE.atomInfo["CD1"].radius = 1.61;
    PHE.atomInfo["CD2"].radius = 1.76;
    PHE.atomInfo["CE1"].radius = 1.76;
    PHE.atomInfo["CE2"].radius = 1.76;
    PHE.atomInfo["CZ"].radius = 1.76;

    auto& MET = types[(AminoAcid)"MET"];
    MET.atomInfo["CB"].radius = 1.88;
    MET.atomInfo["CG"].radius = 1.88;
    MET.atomInfo["SD"].radius = 1.77;
    MET.atomInfo["CE"].radius = 1.88;

    auto& ILE = types[(AminoAcid)"ILE"];
    ILE.atomInfo["CB"].radius = 1.88;
    ILE.atomInfo["CG1"].radius = 1.88;
    ILE.atomInfo["CG2"].radius = 1.88;
    ILE.atomInfo["CD1"].radius = 1.88;

    auto& ASP = types[(AminoAcid)"ASP"];
    ASP.atomInfo["CB"].radius = 1.88;
    ASP.atomInfo["CG"].radius = 1.61;
    ASP.atomInfo["OD1"].radius = 1.46;
    ASP.atomInfo["OD2"].radius = 1.42;

    auto& GLU = types[(AminoAcid)"GLU"];
    GLU.atomInfo["CB"].radius = 1.88;
    GLU.atomInfo["CG"].radius = 1.88;
    GLU.atomInfo["CD"].radius = 1.61;
    GLU.atomInfo["OE1"].radius = 1.46;
    GLU.atomInfo["OE2"].radius = 1.42;

    auto& LYS = types[(AminoAcid)"LYS"];
    LYS.atomInfo["CB"].radius = 1.88;
    LYS.atomInfo["CG"].radius = 1.88;
    LYS.atomInfo["CD"].radius = 1.88;
    LYS.atomInfo["CE"].radius = 1.88;
    LYS.atomInfo["NZ"].radius = 1.64;

    auto& ARG = types[(AminoAcid)"ARG"];
    ARG.atomInfo["CB"].radius = 1.88;
    ARG.atomInfo["CG"].radius = 1.88;
    ARG.atomInfo["CD"].radius = 1.88;
    ARG.atomInfo["NE"].radius = 1.64;
    ARG.atomInfo["CZ"].radius = 1.61;
    ARG.atomInfo["NH1"].radius = 1.64;
    ARG.atomInfo["NH2"].radius = 1.64;

    auto& SER = types[(AminoAcid)"SER"];
    SER.atomInfo["CB"].radius = 1.88;
    SER.atomInfo["OG"].radius = 1.46;

    auto& THR = types[(AminoAcid)"THR"];
    THR.atomInfo["CB"].radius = 1.88;
    THR.atomInfo["OG1"].radius = 1.46;
    THR.atomInfo["CG2"].radius = 1.88;

    auto& TYR = types[(AminoAcid)"TYR"];
    TYR.atomInfo["CB"].radius = 1.88;
    TYR.atomInfo["CG"].radius = 1.61;
    TYR.atomInfo["CD1"].radius = 1.76;
    TYR.atomInfo["CD2"].radius = 1.76;
    TYR.atomInfo["CE1"].radius = 1.76;
    TYR.atomInfo["CE2"].radius = 1.76;
    TYR.atomInfo["CZ"].radius = 1.61;
    TYR.atomInfo["OH"].radius = 1.46;

    auto& HIS = types[(AminoAcid)"HIS"];
    HIS.atomInfo["CB"].radius = 1.88;
    HIS.atomInfo["CG"].radius = 1.61;
    HIS.atomInfo["ND1"].radius = 1.64;
    HIS.atomInfo["CD2"].radius = 1.76;
    HIS.atomInfo["CE1"].radius = 1.76;
    HIS.atomInfo["NE2"].radius = 1.64;

    auto& ASN = types[(AminoAcid)"ASN"];
    ASN.atomInfo["CB"].radius = 1.88;
    ASN.atomInfo["CG"].radius = 1.61;
    ASN.atomInfo["OD1"].radius = 1.42;
    ASN.atomInfo["ND2"].radius = 1.64;

    auto& TRP = types[(AminoAcid)"TRP"];
    TRP.atomInfo["CB"].radius = 1.88;
    TRP.atomInfo["CG"].radius = 1.61;
    TRP.atomInfo["CD1"].radius = 1.76;
    TRP.atomInfo["NE1"].radius = 1.61;
    TRP.atomInfo["CE2"].radius = 1.64;
    TRP.atomInfo["CD2"].radius = 1.61;
    TRP.atomInfo["CE3"].radius = 1.76;
    TRP.atomInfo["CZ3"].radius = 1.76;
    TRP.atomInfo["CZ2"].radius = 1.76;
    TRP.atomInfo["CH2"].radius = 1.76;

    auto& ALA = types[(AminoAcid)"ALA"];
    ALA.atomInfo["CB"].radius = 1.88;

    auto& LEU = types[(AminoAcid)"LEU"];
    LEU.atomInfo["CB"].radius = 1.88;
    LEU.atomInfo["CG"].radius = 1.88;
    LEU.atomInfo["CD1"].radius = 1.88;
    LEU.atomInfo["CD2"].radius = 1.88;

    for (auto& [acid, acidInfo]: types) {
        for (auto& [atomType, atomInfo]: acidInfo.atomInfo) {
            acidInfo.heavyAtoms.insert(atomType);
            atomInfo.radius *= angstrom;
        }
    }

    return types;
}