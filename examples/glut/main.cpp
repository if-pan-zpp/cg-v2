#include <mdk/io/seq/LegacyParser.hpp>
using namespace mdk;

int main() {
    auto seq = seq::LegacyParser().read("data/glut.txt");
    return 0;
}