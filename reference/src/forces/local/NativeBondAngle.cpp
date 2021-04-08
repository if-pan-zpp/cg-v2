#include "forces/local/NativeBondAngle.hpp"
using namespace cg::reference;
using namespace std;

NativeBondAngle::NativeBondAngle(PseudoAtoms const &_pseudoAtoms,
                                 NativeStructure const &_ns):
    pseudoAtoms(_pseudoAtoms),
    ns(_ns) {

    // 'enabled' vector specifies for which i's
    // we should calculate bond angle force
    enabled = vector<unsigned char>(pseudoAtoms.n, 0);

    // TODO: calculate 'enabled' based on NativeStructure
    for (size_t i = 0; i + 2 < pseudoAtoms.n; ++i) enabled[i] = 1;

    nativeTheta = vector<Real>(pseudoAtoms.n, 0.);

    // TODO: calculate native_theta
}

void NativeBondAngle::compute(Reals3 &forces) {
    energy = 0.;

    Eigen::Matrix3d d_theta_d_pos;
    Reals3 const &pos = pseudoAtoms.pos;

    for (size_t i = 0; i < pseudoAtoms.n; ++i) {
        if (enabled[i]) {
            Real3 v0 = pos.col(i + 1) - pos.col(i + 0);
            Real3 v1 = pos.col(i + 2) - pos.col(i + 1);

            Real v0_norm = v0.norm();
            Real v1_norm = v1.norm();
            Real3 v0_x_v1 = v0.cross(v1);

            Real theta = acos(-v0.dot(v1) / (v0_norm * v1_norm));

            Real3 grad0 = v0.cross(v0_x_v1).normalized() / v0_norm;
            Real3 grad2 = v1.cross(v0_x_v1).normalized() / v1_norm;

            d_theta_d_pos.col(0) = grad0;
            d_theta_d_pos.col(1) = -grad0 - grad2;
            d_theta_d_pos.col(2) = grad2;

            theta -= nativeTheta[i];
            energy += k * theta * theta;

            Real d_V_d_theta = 2.0 * k * theta;

            forces.block<3,3>(0,i) -= d_V_d_theta * d_theta_d_pos;
        }
    }
}

void NativeBondAngle::dumpResults(Results &results) {
    results.potEnergy += energy;
}
