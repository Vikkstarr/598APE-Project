#include <iostream>
#include <queue>
#include <vector>
#include <string>
#include <thread>
#include <chrono>
#include <unordered_map>
#include "Hasher.h"

// Default hash table parameters
const size_t DEFAULT_TABLE_SIZE = 1000000;
const size_t DEFAULT_MAX_STEPS = 100;

// Helper function to create test k-mers
std::vector<std::string> generateTestKmers(int count, int kmerLength = 31) {
    std::vector<std::string> kmers;
    const char bases[] = "ACGT";
    
    for (int i = 0; i < count; i++) {
        std::string kmer;
        for (int j = 0; j < kmerLength; j++) {
            kmer += bases[rand() % 4];
        }
        kmers.push_back(kmer);
    }
    return kmers;
}

// Helper to populate queue with blocks
void populateQueue(std::queue<KmerBlock*>& queue, 
                   const std::vector<std::string>& kmers, 
                   int blockSize) {
    for (size_t i = 0; i < kmers.size(); i += blockSize) {
        KmerBlock* block = new KmerBlock();
        for (size_t j = i; j < i + blockSize && j < kmers.size(); j++) {
            block->kmers.push_back(kmers[j]);
        }
        queue.push(block);
    }
}

// Manual count for verification
std::unordered_map<std::string, size_t> manualCount(const std::vector<std::string>& kmers) {
    std::unordered_map<std::string, size_t> counts;
    for (const std::string& kmer : kmers) {
        counts[kmer]++;
    }
    return counts;
}

// Compare two maps
bool compareMaps(const std::unordered_map<std::string, size_t>& map1,
                 const std::unordered_map<std::string, size_t>& map2) {
    if (map1.size() != map2.size()) {
        std::cout << "    Size mismatch: " << map1.size() << " vs " << map2.size() << "\n";
        return false;
    }
    
    for (const auto& entry : map1) {
        auto it = map2.find(entry.first);
        if (it == map2.end() || it->second != entry.second) {
            return false;
        }
    }
    return true;
}

void testHasher(unsigned numThreads, int numKmers, int blockSize) {
    std::cout << "\n=== Test: " << numThreads << " thread(s), " 
              << numKmers << " k-mers, block size " << blockSize << " ===\n";
    
    // Generate test data
    std::vector<std::string> testKmers = generateTestKmers(numKmers);
    
    // Create queue and populate it
    std::queue<KmerBlock*> queue;
    populateQueue(queue, testKmers, blockSize);
    
    std::cout << "  Created " << queue.size() << " blocks\n";
    
    // Create hasher with specified threads
    Hasher hasher(queue, numThreads, DEFAULT_TABLE_SIZE, DEFAULT_MAX_STEPS);
    
    // Start timing
    auto start = std::chrono::high_resolution_clock::now();
    
    // Launch worker threads
    std::vector<std::thread> threads;
    for (unsigned i = 0; i < numThreads; i++) {
        threads.push_back(std::thread(&Hasher::worker, &hasher, i));
    }
    
    // Signal completion
    hasher.signalComplete();
    
    // Wait for all threads
    for (std::thread& t : threads) {
        t.join();
    }
    
    // Stop timing
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Merge results
    hasher.mergeResults();
    const std::unordered_map<std::string, size_t>& results = hasher.getResults();
    
    // Verify correctness
    std::unordered_map<std::string, size_t> expected = manualCount(testKmers);
    bool correct = compareMaps(results, expected);
    
    std::cout << "  Time: " << duration.count() << " ms\n";
    std::cout << "  Unique k-mers found: " << results.size() << "\n";
    std::cout << "  Expected unique k-mers: " << expected.size() << "\n";
    std::cout << "  Results match: " << (correct ? "YES ✓" : "NO ✗") << "\n";
    
    if (!correct) {
        std::cout << "  ERROR: Results don't match!\n";
        // Print first few mismatches for debugging
        int mismatches = 0;
        for (const auto& entry : expected) {
            auto it = results.find(entry.first);
            if (it == results.end()) {
                std::cout << "    Missing: " << entry.first << " (expected " << entry.second << ")\n";
                if (++mismatches >= 5) break;
            } else if (it->second != entry.second) {
                std::cout << "    Wrong count for " << entry.first 
                          << ": got " << it->second << ", expected " << entry.second << "\n";
                if (++mismatches >= 5) break;
            }
        }
    }
}

void testDuplicates() {
    std::cout << "\n=== Test: Duplicate k-mer counting ===\n";
    
    std::queue<KmerBlock*> queue;
    
    // Create blocks with duplicate k-mers
    KmerBlock* block1 = new KmerBlock();
    block1->kmers.push_back("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    block1->kmers.push_back("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    block1->kmers.push_back("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    queue.push(block1);
    
    KmerBlock* block2 = new KmerBlock();
    block2->kmers.push_back("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    block2->kmers.push_back("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    queue.push(block2);
    
    Hasher hasher(queue, 2, DEFAULT_TABLE_SIZE, DEFAULT_MAX_STEPS);
    
    std::vector<std::thread> threads;
    for (unsigned i = 0; i < 2; i++) {
        threads.push_back(std::thread(&Hasher::worker, &hasher, i));
    }
    
    hasher.signalComplete();
    
    for (std::thread& t : threads) {
        t.join();
    }
    
    hasher.mergeResults();
    const std::unordered_map<std::string, size_t>& results = hasher.getResults();
    
    std::cout << "  K-mer counts:\n";
    for (const auto& entry : results) {
        std::cout << "    " << entry.first << ": " << entry.second << "\n";
    }
    
    bool correct = (results.size() == 2 &&
                   results.at("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA") == 3 &&
                   results.at("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT") == 2);
    
    std::cout << "  PASS: " << (correct ? "YES ✓" : "NO ✗") << "\n";
}

void testEmptyQueue() {
    std::cout << "\n=== Test: Empty queue ===\n";
    
    std::queue<KmerBlock*> queue;
    Hasher hasher(queue, 2, DEFAULT_TABLE_SIZE, DEFAULT_MAX_STEPS);
    
    std::vector<std::thread> threads;
    for (unsigned i = 0; i < 2; i++) {
        threads.push_back(std::thread(&Hasher::worker, &hasher, i));
    }
    
    hasher.signalComplete();
    
    for (std::thread& t : threads) {
        t.join();
    }
    
    hasher.mergeResults();
    const std::unordered_map<std::string, size_t>& results = hasher.getResults();
    
    bool correct = (results.size() == 0);
    std::cout << "  Empty result: " << (correct ? "YES ✓" : "NO ✗") << "\n";
}

void speedComparison() {
    std::cout << "\n=== Speed Comparison ===\n";
    const int numKmers = 5000000;
    const int blockSize = 100;
    
    std::vector<unsigned> threadCounts = {1, 2, 4, 8};
    std::vector<long long> times;
    
    for (unsigned numThreads : threadCounts) {
        std::vector<std::string> testKmers = generateTestKmers(numKmers);
        std::queue<KmerBlock*> queue;
        populateQueue(queue, testKmers, blockSize);
        
        Hasher hasher(queue, numThreads, DEFAULT_TABLE_SIZE, DEFAULT_MAX_STEPS);
        
        auto start = std::chrono::high_resolution_clock::now();
        
        std::vector<std::thread> threads;
        for (unsigned i = 0; i < numThreads; i++) {
            threads.push_back(std::thread(&Hasher::worker, &hasher, i));
        }
        
        hasher.signalComplete();
        
        for (std::thread& t : threads) {
            t.join();
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        times.push_back(duration.count());
        std::cout << "  " << numThreads << " thread(s): " << duration.count() << " ms\n";
    }
    
    if (times.size() > 1 && times[0] > 0) {
        std::cout << "  Speedup (1 vs " << threadCounts.back() << " threads): " 
                  << (double)times[0] / times.back() << "x\n";
    }
}

int main() {
    std::cout << "=== Hasher Test Suite ===\n";
    
    // Test 1: Single thread
    testHasher(1, 1000, 50);
    
    // Test 2: Multiple threads with small data
    testHasher(2, 1000, 50);
    testHasher(4, 1000, 50);
    
    // Test 3: Multiple threads with larger data
    testHasher(2, 5000, 100);
    testHasher(4, 5000, 100);
    testHasher(8, 5000, 100);
    
    // Test 4: Different block sizes
    testHasher(4, 2000, 10);   // Small blocks
    testHasher(4, 2000, 500);  // Large blocks
    
    // Test 5: Duplicate counting
    testDuplicates();
    
    // Test 6: Empty queue
    testEmptyQueue();
    
    // Test 7: Speed comparison
    speedComparison();
    
    std::cout << "\n=== All Tests Complete ===\n";
    
    return 0;
}