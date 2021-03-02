#include <iostream>
#include "loaders/SequenceFile.hpp"

int main() {
    auto seqfile = loaders::SequenceFile("data/example3/glut.txt");
    return EXIT_SUCCESS;
}
