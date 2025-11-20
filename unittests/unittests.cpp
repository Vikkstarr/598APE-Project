#include <iostream>
#include <vector> 
#include <string>
#include <cassert>
#include <fstream>
#include "../src/phase1.h"
#include "data_structs.h"

// Test

// Super simple version: just read 1 file into bundles
class FastReader {
    std::string path;
    size_t blockSize;

public:
    FastReader(std::string p, size_t bs = 1 << 20)  // default 1MB chunks
        : path(std::move(p)), blockSize(bs) {}

   std::vector<FastBundle> readFile() {
        std::ifstream in(path);
        if (!in) {
            throw std::runtime_error("Could not open file: " + path);
        }

        std::vector<FastBundle> bundles;
        std::string line;
        std::string seqBuffer;
        seqBuffer.reserve(blockSize * 2);

        while (std::getline(in, line)) {

            if (line.empty() || line[0] == '>') continue;

            // Add DNA to buffer
            seqBuffer += line;

            // When buffer exceeds blockSize: push buffer to bundle 
            while (seqBuffer.size() >= blockSize) {
                FastBundle bundle(blockSize);
                bundle.addBlock(seqBuffer.data(), blockSize);
                bundle.finalize();
                bundles.push_back(std::move(bundle));

                seqBuffer.erase(0, blockSize);
            }
        }

        // Push last  bundle
        if (!seqBuffer.empty()) {
            FastBundle bundle(seqBuffer.size());
            bundle.addBlock(seqBuffer.data(), seqBuffer.size());
            bundle.finalize();
            bundles.push_back(std::move(bundle));
        }

        return bundles;
    }


};


// basic stuff here
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

// now nitty gritty
void testSuperMerToKmers() {
    std::string superMer = "AAGAA";
    int k = 3;
    auto kmers = superMerToKmers(superMer, k);
    std::vector<std::string> expected = {"AAG", "AGA", "GAA"};
    assert(kmers == expected);
    std::cout << "testSuperMerToKmers passed.\n";
}

// Making sure reader works for batching 
void testFastReader_Blocking() {
    // Create a temporary FASTA file
    const std::string filename = "temp_test.fasta";
    std::ofstream out(filename);
    out << ">seq1\n";
    out << "AAGTCCGTA\n";
    out << "GGTAC\n";
    out.close();

    // Batching
    FastReader reader("temp_test.fasta");
    auto bundles = reader.readFile();

    // 3 bundles: AAGTC, CGTAG, GTAC
    std::vector<std::string> expected = {"AAGTC", "CGTAG", "GTAC"};
    assert(bundles.size() == expected.size());
    for (size_t i = 0; i < expected.size(); i++) {
        std::string bundleStr(bundles[i].data.begin(), bundles[i].data.end());
        assert(bundleStr == expected[i]);
    }

    std::cout << "testFastReader_Blocking passed.\n";

    // Clean up temporary file
    std::remove(filename.c_str());
}


int main() {
    testGenerateKmers();
    testComputeMinimizer();
    testComputeAllMinimizers();
    testComputeSuperMers();
    testSuperMerToKmers();
    testFastReader_Blocking();

    std::cout << "All tests passed!\n";
    return 0;
}

