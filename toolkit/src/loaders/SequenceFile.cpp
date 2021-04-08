#include "loaders/SequenceFile.hpp"
#include <fstream>
using namespace cg::toolkit;
using namespace std;

SequenceFile::SequenceFile(filesystem::path const &path) {
    auto file = ifstream(path);
    int nchains;
    file >> nchains;

    for (Index i = 0; i < nchains; ++i) {
        auto &chain = model.chains[i];
        NativeStructure ns;

        /* Load residue codes. */
        int nresidues;
        std::string residue_codes;
        file >> nresidues;
        file >> residue_codes;

        for (auto res: residue_codes) {
            chain.residues.emplace_back(res);
        }

        /* Load contact maps. */
        int ncontactMaps;
        file >> ncontactMaps;

        for (Index j = 0; j < ncontactMaps; ++j) {
            /* Use path relative to the parent dir of the sequence file. */
            string contactMapName;
            file >> contactMapName;
            auto contactMapPath = path.parent_path() / contactMapName;

            auto contactMapIter = mapFiles.find(contactMapPath);
            if (contactMapIter == mapFiles.end()) {
                auto contactMapFile = ifstream(contactMapPath);
                mapFiles[contactMapPath] = ContactMapFile(contactMapFile);
            }

            auto const &contactMap = mapFiles[contactMapPath];
            chain.structuredParts.push_back(contactMap.ns);
        }
    }
}
