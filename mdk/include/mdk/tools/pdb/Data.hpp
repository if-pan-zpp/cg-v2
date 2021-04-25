#pragma once
#include <vector>
#include <mdk/tools/pdb/Record.hpp>
#include <mdk/tools/model/AtomicModel.hpp>

namespace mdk::pdb {
    class Data {
    public:
        std::vector<records::Record> records;

    public:
        Data() = default;
        explicit Data(AtomicModel const& model);
        AtomicModel asModel() const;
    };
}
