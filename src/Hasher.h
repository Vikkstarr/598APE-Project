#ifndef HASHER_H
#define HASHER_H

#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <unordered_map>
#include <string>
#include "QuadraticHashTable.h"
#include "data_structs.h"

class Hasher {
private:
    std::queue<KmerBlock*>& inputQueue;
    std::mutex queueLock;
    std::condition_variable cv;
    
    std::unordered_map<std::string, size_t> globalMap;
    
    // Overflow handling
    std::vector<std::string> overflow;
    std::mutex overflowLock;
    
    unsigned numThreads;
    bool workComplete;

public:
    std::vector<QuadraticHashTable> threadTables;  // Made public for debugging access
    
    Hasher(std::queue<KmerBlock*>& queue, unsigned threads, size_t tableSize, size_t maxSteps);
    
    void worker(unsigned threadId);
    void mergeResults();
    void writeResults(std::string filename);
    void signalComplete();
    
    const std::unordered_map<std::string, size_t>& getResults() const;
};

#endif