#include "data/Topology.hpp"
#include <unordered_set>
using namespace cg;
using namespace std;

static bool operator<(pair<int, int> const &p1, pair<int, int> const &p2) {
    return p1.first < p2.first ||
           (p1.first == p2.first && p1.second < p2.second);
}

class Pairs {
public:
    pair<int, int> excl1, overlap, excl2;

    Pairs(pair<int, int> range1, pair<int, int> range2) {
        if (range1.first > range2.first)
            swap(range1, range2);

        overlap = make_pair(range2.first, range1.second);
        excl1 = make_pair(range1.first, range2.first);
        excl2 = make_pair(range1.second, range2.second);
    }
};

Neighborhood::Neighborhood(Topology *top, Neighborhood::Spec spec) {
    this->top = top;
    this->spec = spec;

    auto minDistSq = spec.cutoff + 2.0 * spec.pad;
    minDistSq *= minDistSq;

    for (auto const& [type1, type2]: spec.pairTypes) {
        int exclusionIx = 0;
        auto range1 = top->system->typeRanges.at(type1);
        auto range2 = top->system->typeRanges.at(type2);
        if (range1.first > range2.first)
            swap(range1, range2);
        auto& [start1, end1] = range1;
        auto& [start2, end2] = range2;

        /* Note: naive implementation */
        for (int ix1 = start1; ix1 < end1; ++ix1) {
            auto chain1 = top->ls->chainId.at(ix1);
            auto pos1 = top->system->pos.col(ix1);

            for (int ix2 = (type1 == type2 ? ix1 : start2); ix2 < end2; ++ix2) {
                auto chain2 = top->ls->chainId.at(ix2);
                auto pos2 = top->system->pos.col(ix2);

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

    startPos = top->system->pos;
    startCell = top->cell;
    maxCorrectDist = spec.cutoff + 2.0 * spec.pad;
}

void Neighborhood::updateMaxCorrectDist() {
    maxCorrectDist = spec.cutoff + 2.0 * spec.pad;

    Real maxDisplacement = 0.0;
    unordered_set<string> types;
    for (auto const& [type1, type2]: spec.pairTypes) {
        types.insert(type1);
        types.insert(type2);
    }
    for (auto type: types) {
        auto [start, end] = top->system->typeRanges.at(type);
        for (int ix = start; ix < end; ++ix) {
            auto curPos = top->system->pos.col(ix);
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

Neighborhood const& Topology::createNeighborhood(const Neighborhood::Spec &spec) {
    return neighborhoods.emplace_back(this, spec);
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
    for (auto& neighborhood: neighborhoods)
        neighborhood.update();
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

