#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <thread>
#include "phase1.h"
#include "Hasher.h"
#include "data_structs.h"

#include <queue>
#include <mutex>
#include <condition_variable>


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


std::vector<std::string> generateKmers(const std::string &seq, int k) {
    std::vector<std::string> kmers;
    if ((int)seq.size() < k) return kmers;  // too short

    for (size_t i = 0; i <= seq.size() - k; i++) {
        kmers.push_back(seq.substr(i, k));
    }
    return kmers;
}

std::string computeMinimizer(const std::string &seq, int m, int k) {
    std::string minimizer = seq.substr(0, m);
    for (size_t i = 1; i < seq.size() - m + 1; i++) {
        std::string current_mmer = seq.substr(i, m);
        if (current_mmer < minimizer) {
            minimizer = current_mmer;
        }
    }
    return minimizer;
}

std::vector<std::string> computeAllMinimizers(const std::string &seq, int m, int k) {
    auto kmers = generateKmers(seq, k);
    std::vector<std::string> minimizers(kmers.size());

    for (size_t i = 0; i < kmers.size(); i++) {
        minimizers[i] = computeMinimizer(kmers[i], m, k);
    }

    return minimizers;
}

std::vector<std::string> computeSuperMers(const std::string &seq, int m, int k) {
    auto kmers = generateKmers(seq, k);
    auto minimizers = computeAllMinimizers(seq, m, k);

    std::vector<std::string> superMers;
    std::string curSuperMer = "";
    std::string curMinimizer = "";

    for(size_t i = 0; i < kmers.size(); i++) {
        if(i == 0) {
            curSuperMer = kmers[i];
            curMinimizer = minimizers[i];
            continue;
        }
        if (minimizers[i] == curMinimizer) {
            curSuperMer += kmers[i].back();
        } else {
            superMers.push_back(curSuperMer);
            curSuperMer = kmers[i];
            curMinimizer = minimizers[i];
        }
    }
    superMers.push_back(curSuperMer);
    return superMers;
}

std::vector<std::string> superMerToKmers(const std::string& superMer, int k) {
    std::vector<std::string> kmers;

    if (superMer.size() < (size_t)k) return kmers;

    kmers.reserve(superMer.size() - k + 1);
    for (size_t i = 0; i <= superMer.size() - k; i++) {
        kmers.push_back(superMer.substr(i, k));
    }
    return kmers;
}

void pushSuperMersToQueue(const std::vector<std::string>& superMers, int k,
                           std::queue<KmerBlock*>& inputQueue,
                           std::mutex& queueLock,
                           std::condition_variable& cv) {
    std::vector<KmerBlock*> localBatch;
    localBatch.reserve(10);

    for(const auto& superMer: superMers) {
        auto kmers = superMerToKmers(superMer, k);
        KmerBlock* block = new KmerBlock(kmers.size());
        for(const auto& kmer: kmers) {
            block->kmers.push_back(kmer);
        }
        localBatch.push_back(block);

        if (localBatch.size() >= 10) {
            {
                std::lock_guard<std::mutex> lock(queueLock);
                for(auto b: localBatch) {
                    inputQueue.push(b);
                }
            }
            cv.notify_all();
            localBatch.clear();
        }
    }
    
    if (!localBatch.empty()) {
        {
            std::lock_guard<std::mutex> lock(queueLock);
            for(auto b: localBatch) {
                inputQueue.push(b);     
            }
        }
        cv.notify_all();    
    }
}

#include <random>

void generateTestFasta(const std::string& filename, size_t length) {
    std::ofstream out(filename);
    const char bases[] = {'A', 'C', 'G', 'T'};

    std::mt19937_64 rng(12345);
    std::uniform_int_distribution<int> dist(0, 3);

    out << ">synthetic\n";
    for (size_t i = 0; i < length; i++) {
        out << bases[dist(rng)];
    }
    out << "\n";
}

int main() {
    const int k = 31;
    const int m = 15;
    const unsigned NUM_THREADS = 4;
    const size_t HASH_TABLE_SIZE = 1000000;
    const size_t MAX_PROBE_STEPS = 100;

    std::cout << "Generating test FASTA...\n";
    generateTestFasta("test_large.fasta", 5'000'000);

    std::cout << "Reading FASTA bundles...\n";
    FastReader reader("test_large.fasta");
    auto bundles = reader.readFile();
    std::cout << "Read " << bundles.size() << " bundles\n";

    std::cout << "Computing super-mers...\n";
    std::vector<std::string> superMers;
    for (auto& b : bundles) {
        auto partial = computeSuperMers(
            std::string(b.data.begin(), b.data.end()),
            m,
            k
        );
        superMers.insert(superMers.end(), partial.begin(), partial.end());
    }
    std::cout << "Total super-mers: " << superMers.size() << "\n";

    std::queue<KmerBlock*> inputQueue;
    std::mutex queueLock;
    std::condition_variable cv;

    std::cout << "Initializing Hasher...\n";
    Hasher hasher(inputQueue, NUM_THREADS, HASH_TABLE_SIZE, MAX_PROBE_STEPS);
    
    std::vector<std::thread> threads;
    
    std::cout << "Launching " << NUM_THREADS << " worker threads...\n";
    for (unsigned i = 0; i < NUM_THREADS; i++) {
        threads.emplace_back([&hasher, i, &queueLock, &cv, &inputQueue]() {
            auto& table = hasher.threadTables[i];
            
            while (true) {
                KmerBlock* block = nullptr;
                {
                    std::unique_lock<std::mutex> lock(queueLock);
                    
                    cv.wait(lock, [&inputQueue, &hasher]() { 
                        return !inputQueue.empty() || hasher.workComplete; 
                    });
                    
                    if (!inputQueue.empty()) {
                        block = inputQueue.front();
                        inputQueue.pop();
                    } else if (hasher.workComplete && inputQueue.empty()) {
                        break;
                    }
                }
                
                if (block) {
                    for (const auto& kmer : block->kmers) {
                        table.insert(kmer);
                    }
                    delete block;
                }
            }
        });
    }

    std::cout << "Pushing super-mers to queue...\n";
    pushSuperMersToQueue(superMers, k, inputQueue, queueLock, cv);
    
    std::cout << "Queue size after batching: " << inputQueue.size() << "\n";

    std::cout << "Signaling completion...\n";
    hasher.signalComplete();
    
    std::cout << "Waiting for threads to finish...\n";
    for (auto& t : threads) {
        t.join();
    }

    std::cout << "Merging results from " << NUM_THREADS << " threads...\n";
    hasher.mergeResults();
    
    const auto& results = hasher.getResults();
    std::cout << "Total unique k-mers: " << results.size() << "\n";

    std::cout << "Writing results to output.txt...\n";
    hasher.writeResults("output.txt");

    std::cout << "\n=== Statistics ===\n";
    std::cout << "Total unique k-mers: " << results.size() << "\n";
    
    size_t maxCount = 0;
    std::string maxKmer;
    for (const auto& [kmer, count] : results) {
        if (count > maxCount) {
            maxCount = count;
            maxKmer = kmer;
        }
    }
    std::cout << "Most frequent k-mer: " << maxKmer << " (count: " << maxCount << ")\n";
    
    std::cout << "\nSample k-mers:\n";
    int shown = 0;
    for (const auto& [kmer, count] : results) {
        if (shown++ >= 5) break;
        std::cout << "  " << kmer << " -> " << count << "\n";
    }

    std::cout << "\nProcessing complete!\n";
    return 0;
}