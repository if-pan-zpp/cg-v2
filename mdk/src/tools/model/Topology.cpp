#include "tools/model/Topology.hpp"
using namespace mdk;

void Topology::setCell(const Eigen::Vector3d &newCell) {
    cell = newCell;
    for (int i = 0; i < 3; ++i) {
        if (cell[i] != 0) cellInv[i] = 1.0/cell[i];
        else cellInv[i] = 0.0;
    }
}
