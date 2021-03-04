#pragma once
#include "utils/AminoAcid.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <vector>

namespace CG {
    /* Note: vectors are column vectors. */
    template<int N, typename T>
    using Vec = Eigen::Matrix<T, N, 1>;

    /* Note: columns are vectors, rows are dimensions. */
    template<int N, typename T>
    using VecList = Eigen::Matrix<T, N, Eigen::Dynamic>;

    /* Numerical entities. */
    using Real = double;
    using RealList = VecList<1, Real>;
    using Real3 = Vec<3, Real>;
    using Real3List = VecList<3, Real>;
    using RealAffine3 = Eigen::Transform<Real, 3, Eigen::Affine>;

    /* Note: in 3.4 version of Eigen, there is better indexing facilities. */
    using Index = unsigned;
    using IndexList = std::vector<Index>;
    using IndexPairList = std::vector<std::pair<Index, Index>>;
}
