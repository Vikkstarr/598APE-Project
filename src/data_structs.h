#ifndef DATA_STRUCTS_H
#define DATA_STRUCTS_H

#include <vector>
#include <string>

// Kmer block structure for batch processing
struct KmerBlock {
    std::vector<std::string> kmers;

    KmerBlock(size_t expectedKmers) {
        kmers.reserve(expectedKmers);
    }

    KmerBlock() = default;
};

// A "bundle" is just a block of raw bytes
struct FastBundle {
    std::vector<char> data;
    bool finalized = false;

    FastBundle(size_t size) { data.reserve(size); }

    void addBlock(const char* buf, size_t n) {
        data.insert(data.end(), buf, buf + n);
    }

    void finalize() { finalized = true; }
};

#endif