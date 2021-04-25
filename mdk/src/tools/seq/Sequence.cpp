#include "tools/seq/Sequence.hpp"
using namespace std;
using namespace mdk;
using namespace mdk::seq;

CoarseModel Sequence::asModel() const {
    CoarseModel model;
    std::unordered_map<std::string, int> cmapIdxMap;

    for (auto const& [cmapName, cmap]: contactMaps) {
        auto& sp = model.addSP();
        cmapIdxMap[cmapName] = sp.idx;
        sp.off = cmap.offset;
        sp.len = cmap.len;
        sp.angle = cmap.angle;
        sp.dihedral = cmap.dihedral;
    }

    for (auto const& chain: chains) {
        auto& chainThere = model.addChain();

        for (auto const& code: chain.codes) {
            auto& resThere = model.addResidue(&chainThere);
            resThere.type = (string)(AminoAcid)code;
            resThere.pos = { 0.0, 0.0, 0.0 };
        }

        int n = chainThere.last - chainThere.first;
        chainThere.tethers = vector<double>(n, 3.8*angstrom);

        for (auto const& cmapName: chain.contactMaps) {
            auto const& sp = contactMaps.at(cmapName);
            auto& spThere = model.structuredParts[cmapIdxMap.at(cmapName)];
            chainThere.structuredParts.push_back(spThere.idx);

            for (auto const& cont: sp.contacts) {
                auto& contThere = model.addContact();
                contThere.type = "CMAP";
                contThere.dist0 = cont.dist0;
                for (int i = 0; i < 2; ++i) {
                    contThere.res[i] = chainThere.first + sp.offset + cont.res[i];
                }
            }
        }
    }

    return model;
}
