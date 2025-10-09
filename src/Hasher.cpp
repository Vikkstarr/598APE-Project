#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <algorithm>

class Hasher {
    private:
        unsigned numThreads;
        bool workComplete = false; // might have to make this an atomic

        std::mutex queueLock;
        std::queue<KmerBlock*>& inputQueue;

        std::vector<std::unordered_map<std::string, size_t>> threadMaps;

    public:

    Hasher(std::queue<KmerBlock*>& queue)
        : inputQueue(queue) {}
    
    void worker(unsigned threadId) {
        auto& map = threadMaps[threadId];
    
        while (true) {
            KmerBlock* block = nullptr;
            std::unique_lock<std::mutex> lock(queueLock);
            
            if (!inputQueue.empty()) {
                block = inputQueue.front();
                inputQueue.pop();
                lock.unlock();
            } else if (workComplete) {
                break;
            }

            if (block) {
                for (const auto& kmer : batch->kmers) {
                    map[kmer]++;
                }
            } else {
                // maybe add cv instead
                std::this_thread::sleep_for(std::chrono::milliseconds(10)); // sleep for a bit because there's no work that can be processed
            }
        }
    }

    void mergeResults() {
        std::unordered_map<std::string, size_t> globalMap;
        for (const auto& threadMap : threadMaps) {
            for (const auto& entry : threadMap) {
                globalMap[entry.first] += entry.second;
            }
        }
    }

}

class QuadraticHashTable {
    private:
        std::vector<std::string> keys;
        std::vector<size_t> values;
        size_t tableSize;
        size_t numElements;
        size_t maxSteps;  // Max probes before giving up (could also do beta (fill ratio))

    public:
        QuadraticHashTable(size_t size = 1009, size_t maxSteps = 5) 
            : tableSize(size), numElements(0), maxSteps(maxSteps) {
                keys.resize(size);
                values.resize(size);

                for (int i = 0; i < tableSize; i++) {
                    keys[i] = "";
                    values[i] = 0;
                }
        }

        bool insert(const std::string* kmer) {
            size_t hashPos = kmer.hash() % tableSize;
            size_t i = 0;
            
            while (true) {
                hashPos = (hashPos + i * i) % tableSize;

                if (keys[hashPos] == "") {
                    keys[hashPos] = kmer;
                    values[hashPos] = 1;
                    numElements++;
                    return true;
                }
                
                if (keys[hashPos] == kmer) {
                    values[hashPos]++;
                    return true;
                }
                
                // collision
                ++i;
                if (i > maxSteps) {
                    return false;
                }
            }
        }

        template<int k> uint64_t
        hash(const std::string& x, const int i) {
            // interpret x as sequence of 32 bases
            uint64_t* data = (uint64_t*) &x;
            uint64_t res = 0; // hash value

            // for each block of 32 bases
            for(int c = 0; c< (k + 31) / 32; ++c) res += 453569 * res + data[c];

            return res + 5696063 * i * i;
        }
}