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
#include "data_structs.h"

class Hasher {
    private:
        std::queue<KmerBlock*>& inputQueue;
        unsigned numThreads;
        
        std::mutex queueLock;
        std::condition_variable cv;
        
        std::unordered_map<std::string, size_t> globalMap;

    public:
        std::vector<QuadraticHashTable> threadTables;
        bool workComplete;

        Hasher(std::queue<KmerBlock*>& queue, unsigned threads, size_t tableSize, size_t maxSteps);
        
        void worker(unsigned threadId);
        void mergeResults();
        void writeResults(std::string filename);
        void signalComplete();
        
        const std::unordered_map<std::string, size_t>& getResults() const;
};

#endif