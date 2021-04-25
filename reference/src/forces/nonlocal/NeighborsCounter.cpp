#include "forces/nonlocal/NeighborsCounter.hpp"
using namespace cg::reference;

#include <iostream>
using namespace std;

NeighborsCounter::NeighborsCounter(PseudoAtoms const &pseudoAtoms,
                                   NativeStructure const &ns,
                                   Topology const &top,
                                   Neighborhood const &verletList,
                                   SharedData &sharedData):
    pseudoAtoms(pseudoAtoms),
    top(top),
    verletList(verletList),
    sharedData(sharedData),
    neiCount(sharedData.neiCount.size(), 0) {


    for (NativeStructure::Contact const &cnt : ns.contacts) {
        nativeContacts.emplace_back(cnt.residues);
    }
}

void NeighborsCounter::checkPair(int i, int j) {
    Real dist_sq = top.offset(pseudoAtoms.pos.col(i),
                              pseudoAtoms.pos.col(j)).squaredNorm();

    if (dist_sq < cutoffSquared) {
        neiCount[i]++;
        neiCount[j]++;
    }
}

void NeighborsCounter::compute(Reals3 &) {
    for (size_t i = 0; i < pseudoAtoms.n; ++i) neiCount[i] = 0;

    for (auto const &[i, j] : verletList.pairs) {
        checkPair(i, j);
    }

    Integers const &chainId = pseudoAtoms.chainId;

    for (size_t i = 0; i + 1 < pseudoAtoms.n; ++i) {
        if (chainId[i] == chainId[i + 1]) {
            checkPair(i, i + 1);
        }
    }

    for (size_t i = 0; i + 2 < pseudoAtoms.n; ++i) {
        if (chainId[i] == chainId[i + 2]) {
            checkPair(i, i + 2);
        }
    }

    for (pair<int, int> const &cnt : nativeContacts) {
        checkPair(cnt.first, cnt.second);
    }
}

void NeighborsCounter::dumpResults(Results &) {
    int tot = 0;
    for (unsigned x : neiCount) tot += x;
    cout << "         " << tot << endl;
    swap(sharedData.neiCount, neiCount);
}
