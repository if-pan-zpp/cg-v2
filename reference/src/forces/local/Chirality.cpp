#pragma once
#include <forces/local/Chirality.hpp>
using namespace cg;

Chirality::Chirality(System *system, NativeStructure *ns, Real echi) {
    this->system = system;
    this->ns = ns;
    this->echi = echi;
}

void Chirality::compute(Reals *energy, Reals3 *force) {
    /* TODO: implement */
}