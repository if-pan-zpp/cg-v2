#include "loaders/SequenceFile.hpp"
#include <fstream>
using namespace CG;
using namespace std;

SequenceFile::SequenceFile(filesystem::path const& path) {
    auto file = ifstream(path);
    int nchains;
    file >> nchains;

    chains = vector<Chain>(nchains);

    for (int i = 0; i < nchains; ++i) {
        auto& chain = chains[i];

        /* Load residue codes. */
        int nresidues;
        file >> nresidues;
        chain.residue_codes.resize(nresidues);
        file >> chain.residue_codes;

        /* Load contact maps. */
        int ncontact_maps;
        file >> ncontact_maps;
        chain.contact_files = vector<ContactFile*>(ncontact_maps);

        for (int j = 0; j < ncontact_maps; ++j) {
            /* Use path relative to the parent dir of the sequence file. */
            string contact_map_name;
            file >> contact_map_name;
            auto contact_map_path = path.parent_path() / contact_map_name;

            auto contact_map_iter = contact_files.find(contact_map_path);
            if (contact_map_iter == contact_files.end()) {
                /* If new contact map, insert into the map. */
                auto contact_map_file = ifstream(contact_map_path);
                auto& contact_file = contact_files[contact_map_path];
                contact_file = ContactFile(contact_map_file);
                chain.contact_files[j] = &contact_file;
            }
            else {
                /* Otherwise, use existing contact map. */
                chain.contact_files[j] = &contact_map_iter->second;
            }
        }
    }
}
