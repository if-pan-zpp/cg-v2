#include "forces/local/ComplexNativeDihedral.hpp"
using namespace cg::reference;

ComplexNativeDihedral::ComplexNativeDihedral(PseudoAtoms const &pseudoAtoms,
                                             NativeStructure const &ns) :
    pseudoAtoms(pseudoAtoms) {

    // 'enabled' vector specifies for which i's
    // we should calculate dihedral angle force
    enabled = vector<unsigned char>(pseudoAtoms.n, 0);
    nativePhi = vector<Real>(pseudoAtoms.n, 0.);

    // TODO: do this in the general case
    for (size_t i = 0; i + 3 < pseudoAtoms.n; ++i) {
        enabled[i] = 1;
        nativePhi[i] = ns.dihedral(i + 2);
    }
}

void ComplexNativeDihedral::compute(Reals3 &forces) {
    Reals3 const &pos = pseudoAtoms.pos;
    energy = 0.0;
    Eigen::Matrix<Real, 3, 4> d_phi_d_pos; 

    for (size_t i = 0; i + 3 < pseudoAtoms.n; ++i) {
        if (enabled[i]) {
            Real3 v0 = pos.col(i + 1) - pos.col(i + 0);
            Real3 v1 = pos.col(i + 2) - pos.col(i + 1);
            Real3 v2 = pos.col(i + 3) - pos.col(i + 2);

            Real3 normal_012 = v0.cross(v1);
            Real3 normal_123 = v1.cross(v2);
            Real normal_012_len_sq = normal_012.squaredNorm();
            Real normal_123_len_sq = normal_123.squaredNorm();

            Real phi = acos(normal_012.dot(normal_123) /
                              sqrt(normal_012_len_sq * normal_123_len_sq));
            if (normal_012.dot(v2) < 0.) phi = -phi;

            Real v1_len_sq = v1.squaredNorm();
            Real v1_len = sqrt(v1_len_sq);

            d_phi_d_pos.col(0) = -normal_012 * v1_len / normal_012_len_sq;
            d_phi_d_pos.col(3) = normal_123 * v1_len / normal_123_len_sq;

            // TODO: figure out why this.
            Real aux1 = -v0.dot(v1);
            Real aux2 = -v2.dot(v1);

            Real3 df = d_phi_d_pos.col(0) * aux1 - d_phi_d_pos.col(3) * aux2;
            df /= v1_len_sq;

            d_phi_d_pos.col(1) = -d_phi_d_pos.col(0) + df;
            d_phi_d_pos.col(2) = -d_phi_d_pos.col(3) - df;

            phi -= nativePhi[i];

            energy += K1 * (1.0 - cos(phi)) + K3 * (1.0 - cos(3. * phi));
            Real d_V_d_phi = K1 * sin(phi) + K3 * 3. * sin(3. * phi);

            forces.block<3,4>(0,i) -= d_V_d_phi * d_phi_d_pos;
        }
    }
}

void ComplexNativeDihedral::dumpResults(Results &results) {
    results.potEnergy += energy;
}
