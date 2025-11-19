// unit tests for phase1.cpp
#include <iostream>
#include <vector> 
#include <string>
#include <cassert>
#include "../src/phase1.h" // include the implementation file to be tested


// basic for testing 
void testGenerateKmers() {
    std::string seq = "AAGTC";
    int k = 3;
    auto kmers = generateKmers(seq, k);
    std::vector<std::string> expected = {"AAG", "AGT", "GTC"};
    assert(kmers == expected);
    std::cout << "testGenerateKmers passed.\n";
}
void testComputeMinimizer() {
    std::string kmer = "AAGTC";
    int m = 3;
    int k = 5;
    auto minimizer = computeMinimizer(kmer, m, k);
    std::string expected = "AAG";
    assert(minimizer == expected);
    std::cout << "testComputeMinimizer passed.\n";
}   

void testComputeAllMinimizers() {
    std::string seq = "AAGTC";
    int m = 3;
    int k = 5;
    auto minimizers = computeAllMinimizers(seq, m, k);
    std::vector<std::string> expected = {"AAG"};
    assert(minimizers == expected);
    std::cout << "testComputeAllMinimizers passed.\n";
}
void testComputeSuperMers() {
    std::string seq = "AAGAACT";
    int m = 3;
    int k = 5;
    auto superMers = computeSuperMers(seq, m, k);
    std::vector<std::string> expected = {"AAGAA", "ACT"};
    assert(superMers == expected);
    std::cout << "testComputeSuperMers passed.\n";
}
int main() {
    testGenerateKmers();
    testComputeMinimizer();
    testComputeAllMinimizers();
    testComputeSuperMers();
    std::cout << "All tests passed!\n";
    return 0;
}
