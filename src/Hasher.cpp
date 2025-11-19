#include "Hasher.h"
#include <algorithm>
#include <chrono>
#include <fstream>

Hasher::Hasher(std::queue<KmerBlock*>& queue, unsigned threads, size_t tableSize, size_t maxSteps)
    : inputQueue(queue), numThreads(threads), workComplete(false) {
    for (unsigned i = 0; i < numThreads; i++) {
        threadTables.push_back(QuadraticHashTable(tableSize, maxSteps));
    }
}

void Hasher::worker(unsigned threadId) {
    QuadraticHashTable& table = threadTables[threadId];

    while (true) {
        KmerBlock* block = nullptr;
        std::unique_lock<std::mutex> lock(queueLock);
        
        cv.wait(lock, [this]() { return !inputQueue.empty() || workComplete; });
        
        if (!inputQueue.empty()) {
            block = inputQueue.front();
            inputQueue.pop();
            lock.unlock();
        } else if (workComplete) {
            break;
        }

        if (block) {
            for (const auto& kmer : block->kmers) {
                table.insert(kmer);
            }
        }
    }
}

void Hasher::mergeResults() {
    globalMap.clear();
    
    for (const auto& table : threadTables) {
        table.exportToMap(globalMap);
    }
}

void Hasher::writeResults(std::string filename) {
    std::ofstream out(filename);
    for (const auto& [kmer, count] : globalMap) {
        out << kmer << '\t' << count << '\n';
    }
}

void Hasher::signalComplete() {
    {
        std::lock_guard<std::mutex> lock(queueLock);
        workComplete = true;
    }
    cv.notify_all();
}

const std::unordered_map<std::string, size_t>& Hasher::getResults() const {
    return globalMap;
}