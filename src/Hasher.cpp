#include "Hasher.h"
#include <fstream>
#include <iostream>

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
        {
            std::unique_lock<std::mutex> lock(queueLock);
            cv.wait(lock, [this]() { 
                return !inputQueue.empty() || workComplete; 
            });
            
            if (!inputQueue.empty()) {
                block = inputQueue.front();
                inputQueue.pop();
            } else if (workComplete && inputQueue.empty()) {
                break;
            }
        }
        
        if (block) {
            for (const auto& kmer : block->kmers) {
                if (!table.insert(kmer)) {
                    // Insertion failed, add to overflow
                    std::lock_guard<std::mutex> lock(overflowLock);
                    overflow.push_back(kmer);
                }
            }
            delete block;
        }
    }
}

void Hasher::mergeResults() {
    globalMap.clear();
    
    for (const auto& table : threadTables) {
        table.exportToMap(globalMap);
    }
    
    std::unordered_map<std::string, size_t> overflowCounts;
    for (const auto& k : overflow) {
        overflowCounts[k]++;
    }
    
    for (const auto& [k, c] : overflowCounts) {
        globalMap[k] += c;
    }
    
    overflow.clear();
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