#pragma once
#include "math/Types.hpp"

namespace math {
    /* This specifier combo should allow the compiler to inline the values
     * at compile time. */
#define Unit inline constexpr Scalard

    /* Convert from one unit to another. */
#define convert(x, unit1, unit2) ((x*unit1)/unit2)

    /* Distance */
    Unit nm = 1.0;
    Unit angstrom = nm / 1.0e1;
    Unit pm = nm / 1.0e3;
    Unit fm = nm / 1.0e6;
    Unit m = nm * 1.0e9;

    /* Time */
    Unit ns = 1.0;
    Unit ps = ns / 1.0e3;
    Unit fs = ns / 1.0e6;
    Unit s = ns * 1.0e9;

    /* Quantity */
    Unit atom = 1.0;
    Unit NA = 1.0;
    Unit mol = 6.02214076e23 / NA;

    /* Energy */
    Unit eps = 1.0; /* \approx 1.5kcal/mol */
    Unit kcal = eps * mol / 1.5;
    Unit J = kcal / 4184.0;
    Unit eV = 1.602176634e-19 * J;

    /* Temperature */
    Unit eps_over_kB = 1.0;
    Unit kB = eps;
    Unit K = 1.380649e-23 * J / kB;

    /* Mass */
    Unit u = 1.0;
    Unit kg = u * mol / 0.99999999965e-3;

    /* EM stuff */
    Unit e = 1.0;
    Unit C = e / 1.602176634e-19;
    Unit A = C / s;
    Unit c = 299792458.0 * m / s;
    Unit H = kg*m*m/(s*s*A*A);
    Unit mu_0 = 1.25663706212e-6 * H / m;
    Unit varepsilon_0 = 1.0 / (mu_0 * c*c);
}