#pragma once
#include <memory>
#include "tools/pdb/Field.hpp"
#include "tools/pdb/Record.hpp"

namespace mdk::pdb {
    using namespace records;

    class RecordParser {
    protected:
        Record record;

    public:
        std::vector<std::shared_ptr<Field>> fields;

        Record tryParse(std::string const& s) const;
        std::ostream& write(std::ostream& os, Record const& other);
    };

    class AtomParser: public RecordParser {
    public:
        AtomParser();
    };

    class HetatmParser: public RecordParser {
    public:
        HetatmParser();
    };

    class SSBondParser: public RecordParser {
    public:
        SSBondParser();
    };

    class Cryst1Parser: public RecordParser {
    public:
        Cryst1Parser();
    };

    class EndParser: public RecordParser {
    public:
        EndParser();
    };

    class TerParser: public RecordParser {
    public:
        TerParser();
    };
}