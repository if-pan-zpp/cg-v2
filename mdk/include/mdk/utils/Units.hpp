#include <math.h>

namespace mdk
{
#define Unit inline constexpr double

    /* Distance */
    Unit nanometer = 1.0;
    Unit angstrom = nanometer / 1.0e1;
    Unit picometer = nanometer / 1.0e3;
    Unit femtometer = nanometer / 1.0e6;
    Unit meter = nanometer * 1.0e9;
    Unit f77unit = 5.0 * angstrom;

    /* Time */
    Unit tau = 1.0;
    Unit nanosecond = 1.0;
    Unit picosecond = nanosecond / 1.0e3;
    Unit femtosecond = nanosecond / 1.0e6;
    Unit second = nanosecond * 1.0e9;

    /* Quantity */
    Unit atom = 1.0;
    Unit NAvogadro = 1.0;
    Unit mol = 6.02214076e23 / NAvogadro;

    /* Energy */
    Unit eps = 1.0; /* \approx 1.5kcal/mol */
    Unit kcal = eps * mol / 1.5;
    Unit Joule = kcal / 4184.0;
    Unit eV = 1.602176634e-19 * Joule;

    /* Temperature */
    Unit epsDivkB = 1.0;
    Unit kB = eps;
    Unit Kelvin = 1.380649e-23 * Joule / kB;

    /* Mass */
    Unit dalton = 1.0;
    Unit kilogram = dalton * mol / 0.99999999965e-3;

    /* EM stuff */
    Unit echarge = 1.0;
    Unit Coulomb = echarge / 1.602176634e-19;
    Unit Ampere = Coulomb / second;
    Unit cspeed = 299792458.0 * meter / second;
    Unit Henry = kilogram * meter * meter / (second * second * Ampere * Ampere);
    Unit mu_0 = 1.25663706212e-6 * Henry / meter;
    Unit varepsilon_0 = 1.0 / (mu_0 * cspeed * cspeed);

    /* Degrees */
    Unit radian = 1.0;
    Unit degree = (2.0 * M_PI / 360.0) * radian;

#undef Unit
}