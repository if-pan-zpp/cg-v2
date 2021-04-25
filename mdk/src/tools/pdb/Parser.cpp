#include "tools/pdb/Parser.hpp"
#include "utils/Text.hpp"
#include "RecordParsers.hpp"
using namespace mdk::pdb::records;
using namespace mdk::pdb;
using namespace std;

Parser::Parser() {
    parsers = {
        { ((Record)Atom()).index(), make_shared<AtomParser>() },
        { ((Record)SSBond()).index(), make_shared<SSBondParser>() },
        { ((Record)Cryst1()).index(), make_shared<Cryst1Parser>() },
        { ((Record)Ter()).index(), make_shared<TerParser>() },
        { ((Record)End()).index(), make_shared<EndParser>() },
    };
}

Data Parser::read(std::istream& is) {
    Data data;

    for (string line; getline(is, line); ) {
        if (line.size() < 80)
            line.insert(line.size(), 80 - line.size(), ' ');

        for (auto const& [idx, parser]: parsers) {
            auto record = parser->tryParse(line);
            if (!holds_alternative<monostate>(record)) {
                data.records.push_back(record);

                if (holds_alternative<End>(record)) return data;
                else break;
            }
        }
    }

    return data;
}

std::ostream &Parser::write(ostream &os, const Data &data) {
    for (auto const& record: data.records) {
        parsers[record.index()]->write(os, record);
        os << endl;
    }
    return os;
}
