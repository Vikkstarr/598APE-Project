#ifndef HASHER_H
#define HASHER_H

#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <unordered_map>
#include <condition_variable>
#include "QuadraticHashTable.h"

struct KmerBlock {
    std::vector<std::string> kmers;
};

class Hasher {
    private:
        unsigned numThreads;
        bool workComplete;

        std::mutex queueLock;
        std::condition_variable cv;
        std::queue<KmerBlock*>& inputQueue;

        std::vector<QuadraticHashTable> threadTables;
        std::unordered_map<std::string, size_t> globalMap;

    public:
        Hasher(std::queue<KmerBlock*>& queue, unsigned threads = std::thread::hardware_concurrency(), 
               size_t tableSize = 1000003, size_t maxSteps = 10);
        
        void worker(unsigned threadId);
        void mergeResults();
        void signalComplete();
        const std::unordered_map<std::string, size_t>& getResults() const;
};

#endif