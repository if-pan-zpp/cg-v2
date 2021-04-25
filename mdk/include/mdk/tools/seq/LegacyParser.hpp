#include <filesystem>
#include <optional>
#include <mdk/utils/AminoAcid.hpp>
#include <mdk/tools/seq/Sequence.hpp>

namespace mdk::seq {
    class LegacyParser {
    public:
        Sequence read(std::filesystem::path const& path);
        void write(std::filesystem::path const& path, Sequence const& seq);
    };
}