#include "model/Chain.hpp"

using namespace cg::toolkit;

void Chain::intoLine(Real3 start, Real3 direction, Real step) {
    direction = direction.normalized() * step;
    for (Index i = 0; i < positions.cols(); ++i) {
        positions.col(i) = start + direction * i;
    }
}

void Chain::intoSAW(bool dense, bool use_pbc, Real initial_density, Real cutoff, RNG &rng) {
    size_t n = residues.size();
    Real x_min, y_min, z_min, x_max, y_max, z_max;
    if(dense) {
        Real start_volume = n / initial_density;
        Real start_box_size = 0.5 * pow(start_volume, 1. / 3.);
        x_min = y_min = z_min = - start_box_size;
        x_max = y_max = z_max = start_box_size;
        if(use_pbc) {
            cell(0) = x_max - x_min;
            cell(1) = y_max - y_min;
            cell(2) = z_max - z_min;
            cell_inv(0) = 1. / cell(0);
            cell_inv(1) = 1. / cell(1);
            cell_inv(2) = 1. / cell(2);
        }
    } else {
        x_min = y_min = z_min = x_max = y_max = z_max = 0.;
    }
    bool conf_correct = false;


    // TODO: we should only do this in the structured case
    bond_length = 0.0;
    for (size_t i = 0; i + 1 < n; ++i) {
        bond_length += (positions.col(i + 1) - positions.col(i)).norm();
    }
    bond_length /= n - 1;

    Reals3_3 T;
    T.resize(n, Real3_3::Zero());
    Real3List R = Real3List::Zero(3, n);
    Real3 r;
    RealList phi = RealList::Zero(n);
    RealList theta = RealList::Zero(n);

    phi[0] = M_PI / 2.;
    theta[0] = 0.;
    phi[1] = 0.;
    theta[1] = rng.uniform() * M_PI / 3.;
    for(size_t i = 2; i + 1 < n; i++) {
        phi[i] = (2. * rng.uniform() - 1.) * M_PI;
        theta[i] = rng.uniform() * M_PI / 3.;
    }

    for(size_t i = 0; i + 1 < n; i++) {
        T[i] << cos(theta[i]), sin(theta[i]), 0., 
                sin(theta[i]) * cos(phi[i]), -cos(theta[i]) * cos(phi[i]), sin(phi[i]),
                sin(theta[i]) * sin(phi[i]), -cos(theta[i]) * sin(phi[i]), -cos(phi[i]);
    }

    Real theta0 = acos(1. - 2. * rng.uniform());
    Real phi0 = 2. * M_PI * rng.uniform();

    for(size_t i = 0; i + 1 < n; i++) {
        Real3 r;
        r << bond_length * sin(theta0) * cos(phi0),
             bond_length * sin(theta0) * sin(phi0),
             bond_length * cos(theta0);
        for(int j = i; j >= 0; j--) {
            R.col(i) = T[j] * r;
            r = R.col(i);
        }
    }

    Real3 ran;
    ran << (x_max - x_min) * rng.uniform() + x_min, 
           (y_max - y_min) * rng.uniform() + y_min,
           (z_max - z_min) * rng.uniform() + z_min;

    for(size_t i = 0; i < n; i++) {
        positions.col(i) = ran;
        for(size_t j = 0; j < i; j++)
            positions.col(i) += R.col(j);
    }

    for(size_t i = 0; i + 3 < n; i++) {
        for(size_t j = i + 3; j < n; j++) {
            Real3 diff_vec = positions.col(j) - positions.col(i);
            if(use_pbc) {
                for (int d = 0; d < 3; ++d) 
                    diff_vec(d) -= cell(d)*(int)(diff_vec(d)*cell_inv(d));
            }
            Real sq_dist = diff_vec.squaredNorm();
            if(sq_dist < cutoff * cutoff) {
                SAW_attempts_left--;
                if(SAW_attempts_left > 0) { // creation failed, retry
                    intoSAW(dense, use_pbc, initial_density, cutoff, rng);
                    return;
                }
                // TODO announce failure
                return;
            }
        }
    }

    if(!use_pbc) {
        for(size_t i = 0; i < n; i++) {
            x_min = min(x_min, positions.col(i)[0]);
            x_max = max(x_max, positions.col(i)[0]);
            y_min = min(y_min, positions.col(i)[1]);
            y_max = max(y_max, positions.col(i)[1]);
            z_min = min(z_min, positions.col(i)[2]);
            z_max = max(z_max, positions.col(i)[2]);
        }
    }

}

void Chain::apply(RealAffine3 const &aff) {
    for (Index i = 0; i < positions.cols(); ++i) {
        positions.col(i) = aff * positions.col(i);
    }
}
