#include "tools/pdb/Data.hpp"
#include <unordered_map>
#include <unordered_set>
#include "utils/TupleHash.hpp"
using namespace mdk::pdb::records;
using namespace mdk::pdb;
using namespace mdk;
using namespace std;

AtomicModel Data::asModel() const {
    AtomicModel model;
    unordered_set<char> terminated;

    for (auto const& record: records) {
        if (holds_alternative<Atom>(record)) {
            auto& data = get<Atom>(record);

            if (model.chains.count(data.chainID) == 0) {
                model.addChain(data.chainID);
            }
            auto& chain = model.chains[data.chainID];
            if (terminated.count(data.chainID) > 0) continue;

            if (model.residues.count(data.residueSeqNum) == 0) {
                auto& res = model.addResidue(data.residueSeqNum, &chain);
                res.type = data.residueName;
            }
            auto& res = model.residues[data.residueSeqNum];

            if (model.atoms.count(data.serialNum) == 0) {
                auto& atom = model.addAtom(data.serialNum, &res);
                atom.type = data.atomName;
                atom.pos = data.pos;
            }
        }
        else if (holds_alternative<Ter>(record)) {
            auto& data = get<Ter>(record);
            auto& res = model.residues[data.residueSeqNum];
            assert(data.residueName == res.type);
            terminated.insert(data.chainId);
        }
        else if (holds_alternative<SSBond>(record)) {
            auto& data = get<SSBond>(record);

            auto& ssbond = model.addContact();
            ssbond.type = "SSBOND";
            ssbond.dist0 = data.dist0;

            for (int i = 0; i < 2; ++i) {
                auto& info = data.res[i];
                auto& res = model.residues[info.residueSeqNum];
                ssbond.atom[i] = res.find("SG");
            }
        }
        else if (holds_alternative<Cryst1>(record)) {
            auto& data = get<Cryst1>(record);
            model.cell = data.size;
        }
        else if (holds_alternative<End>(record)) {
            break;
        }
    }

    return model;
}

Data::Data(const AtomicModel &model) {
    if (model.cell) {
        auto record = records.emplace_back(Cryst1());
        auto& data = get<Cryst1>(record);
        data.size = model.cell.value();
        data.angles = { 0.0, 0.0, 0.0 };
    }

    for (auto const& [chainIdx, chain]: model.chains) {
        AtomicModel::Atom *finalAtom;

        for (auto const& res: chain.residues) {
            for (auto const& atom: res->atoms) {
                auto& record = records.emplace_back(Atom());
                auto& data = get<Atom>(record);
                data.atomName = atom->type;
                data.residueName = atom->res->type;
                data.residueSeqNum = atom->res->serial;
                data.serialNum = atom->serial;
                data.pos = atom->pos;
                data.chainID = atom->res->chain->serial;
                data.element = atom->type[0];

                finalAtom = atom;
            }
        }

        auto& record = records.emplace_back(Ter());
        auto& data = get<Ter>(record);
        auto const& finalRes = finalAtom->res;
        data.serialNum = finalAtom->serial + 1;
        data.residueSeqNum = finalRes->serial;
        data.residueName = finalRes->type;
        data.chainId = chain.serial;
    }

    for (auto const& cont: model.contacts) {
        if (cont.type != "SSBOND") continue;

        auto& record = records.emplace_back(SSBond());
        auto& data = get<SSBond>(record);
        data.dist0 = cont.dist0;
        data.serialNum = cont.idx;
        for (int i = 0; i < 2; ++i) {
            data.res[i].chainId = cont.atom[i]->res->chain->serial;
            data.res[i].residueSeqNum = cont.atom[i]->res->serial;
        }
    }

    records.emplace_back(End());
}
