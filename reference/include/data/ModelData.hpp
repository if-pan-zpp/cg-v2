#pragma once
#include "data/PseudoAtoms.hpp"
#include "data/NativeStructure.hpp"
#include "model/Model.hpp"

namespace cg::reference {
    class ModelData {
    public:
        PseudoAtoms pseudoAtoms;
        NativeStructure ns;

        ModelData(toolkit::Model const& model);
    };
}