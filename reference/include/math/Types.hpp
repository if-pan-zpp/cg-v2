#pragma once
#include <Eigen/Dense>

namespace math {
    template<int N, typename T>
    using Vec = Eigen::Matrix<T, N, 1>;

    template<int N, typename T>
    using VecList = Eigen::Array<T, N, Eigen::Dynamic>;

    /* Numerical entities. */
    using Real = double;
    using RealList = VecList<1, Real>;
    using Real3 = Vec<3, Real>;
    using Real3List = VecList<3, Real>;

    /* Helpers - vectors broadcasts individual elements, whereas dims
     * broadcasts entire axes. */
#define vectors(list) (list.colwise().matrix())
#define dims(list) (list.rowwise().matrix())
}
