#include "loaders/SequenceFile.hpp"
#include <fstream>
using namespace CG;
using namespace std;

SequenceFile::SequenceFile(filesystem::path const& path) {
    auto file = ifstream(path);
    int nchains;
    file >> nchains;

    for (Index i = 0; i < nchains; ++i) {
        auto& chain = model.chains[i];
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
        int ncontact_maps;
        file >> ncontact_maps;

        for (Index j = 0; j < ncontact_maps; ++j) {
            /* Use path relative to the parent dir of the sequence file. */
            string contact_map_name;
            file >> contact_map_name;
            auto contact_map_path = path.parent_path() / contact_map_name;

            auto contact_map_iter = map_files.find(contact_map_path);
            if (contact_map_iter == map_files.end()) {
                auto contact_map_file = ifstream(contact_map_path);
                map_files[contact_map_path] = ContactMapFile(contact_map_file);
            }

            auto const& contact_map = map_files[contact_map_path];
            chain.structured_parts.push_back(contact_map.ns);
        }
    }
}
