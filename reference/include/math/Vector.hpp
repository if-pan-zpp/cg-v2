#pragma once
#include <eigen3/Eigen/Core>
#include <vector>

namespace CG {
    /* Note: vectors are column vectors. */
    template<int N, typename T>
    using Vec = Eigen::Matrix<T, N, 1>;

    /* Note: columns are vectors, rows are dimensions. */
    template<int N, typename T>
    using VecList = Eigen::Array<T, N, Eigen::Dynamic>;

    /* Numerical entities. */
    using Real = double;
    using RealList = VecList<1, Real>;
    using Real3 = Vec<3, Real>;
    using Real3List = VecList<3, Real>;

    /* Note: in 3.4 version of Eigen, there is better indexing facilities. */
    using IndexList = std::vector<int>;
    using IndexPairList = std::vector<std::pair<int, int>>;

    /* Helpers - vectors broadcasts individual elements, whereas dims
     * broadcasts entire axes. */
#define vectors(list) (list.colwise().matrix())
#define dims(list) (list.rowwise().matrix())
}
