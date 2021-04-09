#include "data/Topology.hpp"
#include <unordered_set>
using namespace cg::reference;
using namespace std;

static bool operator<(pair<int, int> const &p1, pair<int, int> const &p2) {
    return p1.first < p2.first ||
           (p1.first == p2.first && p1.second < p2.second);
}

Neighborhood::Neighborhood(Topology *top, Neighborhood::Spec spec) {
    this->top = top;
    this->spec = spec;

    auto minDistSq = spec.cutoff + 2.0 * spec.pad;
    minDistSq *= minDistSq;

    for (auto const &[type1, type2]: spec.pairTypes) {
        int exclusionIx = 0;
        auto range1 = top->pseudoAtoms->typeRanges.at(type1);
        auto range2 = top->pseudoAtoms->typeRanges.at(type2);
        if (range1.first > range2.first)
            swap(range1, range2);
        auto &[start1, end1] = range1;
        auto &[start2, end2] = range2;

        /* Note: naive implementation */
        for (int ix1 = start1; ix1 < end1; ++ix1) {
            auto chain1 = top->ns->chainId.at(ix1);
            auto pos1 = top->pseudoAtoms->pos.col(ix1);

            for (int ix2 = (type1 == type2 ? ix1 : start2); ix2 < end2; ++ix2) {
                auto chain2 = top->ns->chainId.at(ix2);
                auto pos2 = top->pseudoAtoms->pos.col(ix2);

                /* Distance exclusion */
                auto dist = top->offset(pos1, pos2);
                if (dist.squaredNorm() > minDistSq)
                    continue;

                /* "Local" exclusion */
                if (chain1 == chain2) {
                    auto diff = ix2 - ix1;
                    if (diff < spec.minLocalOffset || (diff == 4 && !spec.include4))
                        continue;
                }

                /* List exclusion */
                auto cur = make_pair(ix1, ix2);
                while (exclusionIx < spec.exclusions->size() &&
                       spec.exclusions->at(exclusionIx) < cur)
                    ++exclusionIx;
                if (exclusionIx < spec.exclusions->size() &&
                    spec.exclusions->at(exclusionIx) == cur)
                    continue;

                /* Now, it should be in the list */
                pairs.emplace_back(ix1, ix2);
            }
        }
    }

    startPos = top->pseudoAtoms->pos;
    startCell = top->cell;
    maxCorrectDist = spec.cutoff + 2.0 * spec.pad;
}

void Neighborhood::updateMaxCorrectDist() {
    maxCorrectDist = spec.cutoff + 2.0 * spec.pad;

    Real maxDisplacement = 0.0;
    unordered_set<string> types;
    for (auto const &[type1, type2]: spec.pairTypes) {
        types.insert(type1);
        types.insert(type2);
    }
    for (auto type: types) {
        auto [start, end] = top->pseudoAtoms->typeRanges.at(type);
        for (int ix = start; ix < end; ++ix) {
            auto curPos = top->pseudoAtoms->pos.col(ix);
            Real displacement = top->offset(startPos.col(ix), curPos).norm();
            if (maxDisplacement < displacement)
                maxDisplacement = displacement;
        }
    }
    maxCorrectDist -= 2.0 * maxDisplacement;

    Real boxChange = (startCell - top->cell).norm();
    maxCorrectDist -= boxChange;

    if (maxCorrectDist < 0.0)
        maxCorrectDist = 0.0;
}

void Neighborhood::update() {
    updateMaxCorrectDist();
    if (maxCorrectDist < spec.cutoff)
        forceUpdate();
}

void Neighborhood::forceUpdate() {
    *this = Neighborhood(top, spec);
}

Neighborhood const &Topology::createNeighborhood(const Neighborhood::Spec &spec) {
    unique_ptr<Neighborhood> new_nei = make_unique<Neighborhood>(this, spec);
    Neighborhood const &ref = *new_nei;
    neighborhoods.push_back(move(new_nei));
    return ref;
}

Real3 Topology::offset(const Real3 &p, const Real3 &q) const {
    Real3 off = q - p;
    for (int d = 0; d < 3; ++d) {
        if (cell(d) != 0.0) {
            off(d) -= cell(d)*(int)(off(d)*cell_inv(d));
        }
    }
    return off;
}

void Topology::update() {
    for (auto &neighborhood: neighborhoods)
        neighborhood -> update();
}

void Topology::forceUpdate() {
    // TODO
}

void Topology::updateCellSize(Real3 newCell) {
    bool pbcness_changed = false;

    for (int d = 0; d < 3; ++d) {
        if ((cell(d) == 0.0) ^ (newCell(d) == 0.0))
            pbcness_changed = true;

        cell(d) = newCell(d);
        if (cell(d) != 0.0)
            cell_inv(d) = 1.0/cell(d);
    }

    /* If an axis (un)became periodic, we force an update.
     * Note: is it a good idea? */
    if (pbcness_changed)
        forceUpdate();
}

Topology::Topology(PseudoAtoms const &pseudoAtoms, NativeStructure const &ns) {
    this->pseudoAtoms = &pseudoAtoms;
    this->ns = &ns;
}

