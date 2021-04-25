#pragma once
#include <Eigen/Core>

namespace mdk {
    using Integers = Eigen::MatrixXi;
    using Scalars = Eigen::MatrixXd;
    using Vector = Eigen::Vector3d;
    using Vectors = Eigen::MatrixX3d;

#define vec(vs, i) ((vs).row(i).transpose())
}
