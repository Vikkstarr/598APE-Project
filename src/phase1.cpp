#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "phase1.h"

#include <queue>
#include <mutex>
#include <condition_variable>



// phase 2 first part converting to kmer blocks

struct KmerBlock {
    std::vector<std::string> kmers;

    // block for given expected kmers
    KmerBlock(size_t expectedKmers) {
        kmers.reserve(expectedKmers);
    }

    // Default block (empty, no reserve)
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

    //#pragma omp parallel for
    for (size_t i = 0; i <= seq.size() - k; i++) {
        kmers.push_back(seq.substr(i, k));
    }
    return kmers;
}

// "super-mer of a genome sequence is defined as a substring of maximal 
// length such that all k-mers on that substring share the same minimizer. 
// Hereby, a minimizer of a k-mer is defined as a substring of fixed length 
// m that is minimal with respect to some total ordering on strings of length 
// m. Thus, contiguous k-mers of a genome read are joined to a super-mer if they 
// share the same minimizer.”


std::string computeMinimizer(const std::string &seq, int m, int k) {
    // so basically take the first m-mer of the first k-mer and itetate thorugh every one

    

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

    // Parallelize this loop
    #pragma omp parallel for
    for (size_t i = 0; i < kmers.size(); i++) {
        minimizers[i] = computeMinimizer(kmers[i], m, k);
    }

    return minimizers;
}

std::vector<std::string> computeSuperMers(const std::string &seq, int m, int k) {
    auto kmers = generateKmers(seq, k);

    // Time this part
    //double start = omp_get_wtime();

    auto minimizers = computeAllMinimizers(seq, m, k);

    //double end = omp_get_wtime();
    //std::cout << "Time to compute all minimizers: " << (end - start) << " seconds\n";

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


// Assume KmerBlock and FastReader are defined as before

// Push a single super-mer's kmers into the queue
void pushSuperMerToQueue(const std::string& superMer, int k,
                         std::queue<KmerBlock*>& inputQueue,
                         std::mutex& queueLock,
                         std::condition_variable& cv) {
    if (superMer.size() < (size_t)k) return;

    size_t numKmers = superMer.size() - k + 1;
    KmerBlock* block = new KmerBlock(numKmers);
    block->kmers.reserve(numKmers);
    for (size_t i = 0; i < numKmers; ++i)
        block->kmers.push_back(superMer.substr(i, k));

    {
        std::lock_guard<std::mutex> lock(queueLock);
        inputQueue.push(block);
    }
    cv.notify_all();
}

 #include <fstream>
#include <random>
#include <thread>

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
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <atomic>      


int main() {
    const int k = 31;
    const int m = 15;
    const size_t MAX_QUEUE_SIZE = 1000; // flush threshold
    using clock = std::chrono::high_resolution_clock;

    auto total_start = clock::now();

    std::cout << "Reading FASTA bundles...\n";
    FastReader reader("GCA_000002315.5_GRCg6a_genomic.fna");
    auto bundles = reader.readFile();
    std::cout << "Read " << bundles.size() << " bundles\n";

    auto supermer_start = clock::now();

    std::queue<KmerBlock*> inputQueue;
    std::mutex queueLock;
    std::condition_variable cv;

    const unsigned NUM_THREADS = 8;
    std::atomic<size_t> bundleIndex(0);

    // Worker threads
    std::vector<std::thread> workers;
    for (unsigned t = 0; t < NUM_THREADS; t++) {
        workers.emplace_back([&]() {
            while (true) {
                size_t idx = bundleIndex.fetch_add(1);
                if (idx >= bundles.size()) break;

                auto& b = bundles[idx];
                auto superMers = computeSuperMers(
                    std::string(b.data.begin(), b.data.end()), m, k
                );

                for (auto& sm : superMers) {
                    pushSuperMerToQueue(sm, k, inputQueue, queueLock, cv);

                    // Flush queue if it exceeds MAX_QUEUE_SIZE
                    bool flush = false;
                    {
                        std::lock_guard<std::mutex> lock(queueLock);
                        if (inputQueue.size() >= MAX_QUEUE_SIZE) flush = true;
                    }

                    if (flush) {
                        std::lock_guard<std::mutex> lock(queueLock);
                        // std::cout << "Yay! Queue reached " << inputQueue.size()
                        //           << " blocks. Clearing...\n";

                        while (!inputQueue.empty()) {
                            delete inputQueue.front();
                            inputQueue.pop();
                        }
                    }
                }
            }
        });
    }

    // Join threads
    for (auto& th : workers) th.join();

    // Final flush if anything remains
    if (!inputQueue.empty()) {
        std::cout << "Final flush: clearing remaining " << inputQueue.size() << " blocks\n";
        while (!inputQueue.empty()) {
            delete inputQueue.front();
            inputQueue.pop();
        }
    }

    auto supermer_end = clock::now();
    std::cout << "Super-mer + queue batching time: "
              << std::chrono::duration<double>(supermer_end - supermer_start).count()
              << " sec\n";

    auto total_end = clock::now();
    std::cout << "TOTAL runtime: "
              << std::chrono::duration<double>(total_end - total_start).count()
              << " sec\n";

    std::cout << "Queue processing complete.\n";
    return 0;
}



//________________



// // lowkey same as readFile but pushing to queue in batches of 10
// // with locking and cv notifying
// void pushSuperMersToQueue(const std::vector<std::string>& superMers, int k,
//                            std::queue<KmerBlock*>& inputQueue,
//                            std::mutex& queueLock,
//                            std::condition_variable& cv,
//                            unsigned NUM_THREADS = 8) {
//     size_t N = superMers.size();
//     size_t CHUNK = (N + NUM_THREADS - 1) / NUM_THREADS;

//     // Each thread will have its local vector of KmerBlocks
//     std::vector<std::vector<KmerBlock*>> localBlocks(NUM_THREADS);

//     // Spawn threads
//     std::vector<std::thread> workers;
//     workers.reserve(NUM_THREADS);

//     for (unsigned thr = 0; thr < NUM_THREADS; ++thr) {
//         workers.emplace_back([&, thr]() {
//             size_t startIndex = thr * CHUNK;
//             size_t endIndex = std::min(startIndex + CHUNK, N);

//             auto& myOutput = localBlocks[thr];

//             // Rough reserve
//             size_t numSuperMers = (endIndex > startIndex) ? (endIndex - startIndex) : 1;
//             myOutput.reserve(numSuperMers / 4 + 4);

//             // Process assigned super-mers
//             for (size_t i = startIndex; i < endIndex; ++i) {
//                 const std::string& sm = superMers[i];
//                 if (sm.size() < k) continue;

//                 size_t numKmers = sm.size() - k + 1;
//                 KmerBlock* block = new KmerBlock(numKmers);
//                 block->kmers.reserve(numKmers);

//                 for (size_t j = 0; j < numKmers; ++j) {
//                     block->kmers.push_back(sm.substr(j, k));
//                 }

//                 myOutput.push_back(block);
//             }
//         });
//     }

//     // Wait for threads to finish
//     for (auto& t : workers) t.join();

//     // Merge local blocks into the shared queue
//     {
//         std::lock_guard<std::mutex> lock(queueLock);
//         for (auto& vec : localBlocks) {
//             for (KmerBlock* blk : vec) {
//                 inputQueue.push(blk);
//             }
//         }
//     }

//     // Notify any waiting threads
//     cv.notify_all();
// }

// Sequential version: push super-mers to the shared queue safely
// void pushSuperMersToQueue(const std::vector<std::string>& superMers, int k,
//                            std::queue<KmerBlock*>& inputQueue,
//                            std::mutex& queueLock,
//                            std::condition_variable& cv) {
//     for (const auto& sm : superMers) {
//         if (sm.size() < (size_t)k) continue;

//         size_t numKmers = sm.size() - k + 1;
//         KmerBlock* block = new KmerBlock(numKmers);
//         block->kmers.reserve(numKmers);

//         for (size_t i = 0; i < numKmers; ++i) {
//             block->kmers.push_back(sm.substr(i, k));
//         }

//         // Push safely into the shared queue
//         {
//             std::lock_guard<std::mutex> lock(queueLock);
//             inputQueue.push(block);
//         }
//     }

//     // Notify any waiting threads
//     cv.notify_all();
// }
