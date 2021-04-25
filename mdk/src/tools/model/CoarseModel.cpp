#include "tools/model/CoarseModel.hpp"
#include <Eigen/Geometry>
using namespace std;
using namespace mdk;

CoarseModel::Residue& CoarseModel::addResidue(Chain *chain) {
    auto& res = residues.emplace_back();
    res.idx = residues.size()-1;
    if (chain) {
        res.chainIdx = chain->idx;
        if (chain->first < 0) chain->first = 0;
        chain->last = res.idx+1;
    }
    else {
        res.chainIdx = -1;
    }
    return res;
}

CoarseModel::Chain& CoarseModel::addChain() {
    auto& chain = chains.emplace_back();
    chain.idx = chains.size()-1;
    chain.first = chain.last = -1;
    return chain;
}

CoarseModel::Contact &CoarseModel::addContact() {
    auto& cont = contacts.emplace_back();
    cont.idx = contacts.size()-1;
    return cont;
}

CoarseModel::StructuredPart& CoarseModel::addSP() {
    auto& sp = structuredParts.emplace_back();
    sp.idx = structuredParts.size()-1;
    return sp;
}

void CoarseModel::morphIntoLine() {
    double tether0 = 3.8*angstrom;

    Eigen::Vector3d dir = { 0.0, 0.0, tether0 };
    for (auto& res: residues) {
        res.pos = res.idx * dir;
    }
}

static Eigen::Vector3d sampleBox(Eigen::AlignedBox3d const& box, Random& rand) {
    Eigen::Vector3d v;
    for (int i = 0; i < 3; ++i) {
        v[i] = rand.uniform(box.min()[i], box.max()[i]);
    }
    return v;
}

static Eigen::Vector3d sampleSphere(Random& rand) {
    Eigen::Vector3d v;
    for (int i = 0; i < 3; ++i) {
        v[i] = rand.normal();
    }
    return v.normalized();
}

static Eigen::Vector3d findPerp(Eigen::Vector3d const& v) {
    auto perp1 = v.cross(Eigen::Vector3d::UnitX());
    auto perp2 = v.cross(Eigen::Vector3d::UnitY());
    auto& perp = perp1.norm() > perp2.norm() ? perp1 : perp2;
    return perp.normalized();
}

std::optional<Topology> CoarseModel::morphIntoSAW(Random& rand, bool withPBC, double density, double intersectionDist) {
    Eigen::Vector3d minCorner, maxCorner;
    Eigen::AlignedBox3d box;
    Topology top;
    double tether0 = 3.8*angstrom;
    auto pairs = nonlocalPairs();
    auto state0 = residues;

    if (density > 0) {
        auto vol = residues.size() / (density * atom);
        auto a = 0.5 * pow(vol, 1.0/3.0);
        minCorner = -a * Eigen::Vector3d { 1.0, 1.0, 1.0 };
        maxCorner = -minCorner;
    }
    else {
        minCorner = maxCorner = Eigen::Vector3d { 0.0, 0.0, 0.0 };
    }

    box = Eigen::AlignedBox3d(minCorner, maxCorner);
    top.setCell(maxCorner - minCorner);

    for (int attempt = 0; attempt < 9000; ++attempt) {
        bool failed = false;

        for (auto const& chain: chains) {
            auto pos = sampleBox(box, rand);
            auto dir = sampleSphere(rand);

            for (int i = chain.first; i < chain.last; ++i) {
                residues[i].pos = pos;
                pos += tether0 * dir;

                double spread = rand.uniform(M_PI/12.0, M_PI/2.0);
                auto spreadRot = Eigen::AngleAxisd(spread, findPerp(dir));

                double around = rand.uniform(-M_PI, M_PI);
                auto aroundRot = Eigen::AngleAxisd(around, dir);

                dir = aroundRot * spreadRot * dir;
                dir.normalize();
            }
        }

        for (auto& res: residues) {
            if (res.chainIdx < 0)
                res.pos = sampleBox(box, rand);
        }

        for (auto const& [res1, res2]: pairs) {
            Eigen::Vector3d dx;
            if (!withPBC) dx = res1->pos - res2->pos;
            else dx = top.pbc(res1->pos - res2->pos);

            if (dx.norm() < intersectionDist) {
                failed = true;
                break;
            }
        }

        if (!failed)
            return top;
    }

    residues = state0;
    return nullopt;
}

vector<pair<CoarseModel::Residue*, CoarseModel::Residue*>> CoarseModel::nonlocalPairs() {
    vector<pair<Residue*, Residue*>> pairs;
    for (auto& res1: residues) {
        for (auto& res2: residues) {
            if (res1.chainIdx == res2.chainIdx && res1.chainIdx >= 0)
                continue;
            if (abs(res1.idx - res2.idx) < 3)
                continue;

            pairs.emplace_back(&res1, &res2);
        }
    }

    return pairs;
}
