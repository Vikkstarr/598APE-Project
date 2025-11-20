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
# include <atomic>


// struct KmerBlock {
//     std::vector<std::string> kmers;

//     // block for given expected kmers
//     KmerBlock(size_t expectedKmers) {
//         kmers.reserve(expectedKmers);
//     }

//     // Default block (empty, no reserve)
//     KmerBlock() = default;
// };


// // A "bundle" is just a block of raw bytes
// struct FastBundle {
//     std::vector<char> data;
//     bool finalized = false;

//     FastBundle(size_t size) { data.reserve(size); }

//     void addBlock(const char* buf, size_t n) {
//         data.insert(data.end(), buf, buf + n);
//     }

//     void finalize() { finalized = true; }
// };

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

// int main() {
//     const int k = 6;
//     const int m = 5;
//     const unsigned NUM_THREADS = 8;
//     const size_t HASH_TABLE_SIZE = 10000000;
//     const size_t MAX_PROBE_STEPS = 100;

//     std::cout << "Generating test FASTA...\n";
//     generateTestFasta("test_large.fasta", 5'000'000);

//     std::cout << "Reading FASTA bundles...\n";
//     FastReader reader("test_large.fasta");
//     auto bundles = reader.readFile();
//     std::cout << "Read " << bundles.size() << " bundles\n";

//     std::cout << "Computing super-mers...\n";
//     std::vector<std::string> superMers;
//     for (auto& b : bundles) {
//         auto partial = computeSuperMers(
//             std::string(b.data.begin(), b.data.end()),
//             m,
//             k
//         );
//         superMers.insert(superMers.end(), partial.begin(), partial.end());
//     }
//     std::cout << "Total super-mers: " << superMers.size() << "\n";

//     std::queue<KmerBlock*> inputQueue;
//     std::mutex queueLock;
//     std::condition_variable cv;

//     std::cout << "Initializing Hasher...\n";
//     Hasher hasher(inputQueue, NUM_THREADS, HASH_TABLE_SIZE, MAX_PROBE_STEPS);
    
//     std::vector<std::thread> threads;
    
//     std::cout << "Launching " << NUM_THREADS << " worker threads...\n";
//     for (unsigned i = 0; i < NUM_THREADS; i++) {
//         threads.emplace_back(&Hasher::worker, &hasher, i);
//     }

//     std::cout << "Pushing super-mers to queue...\n";
//     pushSuperMersToQueue(superMers, k, inputQueue, queueLock, cv);
    
//     std::cout << "Queue size after batching: " << inputQueue.size() << "\n";

//     std::cout << "Signaling completion...\n";
//     hasher.signalComplete();
    
//     std::cout << "Waiting for threads to finish...\n";
//     for (auto& t : threads) {
//         t.join();
//     }

//     std::cout << "Merging results from " << NUM_THREADS << " threads...\n";

//     // Debug: Check each thread's table before merging
//     for (unsigned i = 0; i < NUM_THREADS; i++) {
//         std::cout << "Thread " << i << " table:\n";
//         hasher.threadTables[i].printStats();
//     }

//     hasher.mergeResults();
    
//     const auto& results = hasher.getResults();
//     std::cout << "Total unique k-mers: " << results.size() << "\n";

//     std::cout << "Writing results to output.txt...\n";
//     hasher.writeResults("output.txt");

//     std::cout << "\n=== Statistics ===\n";
//     std::cout << "Total unique k-mers: " << results.size() << "\n";
    
//     size_t maxCount = 0;
//     std::string maxKmer;
//     for (const auto& [kmer, count] : results) {
//         if (count > maxCount) {
//             maxCount = count;
//             maxKmer = kmer;
//         }
//     }
//     std::cout << "Most frequent k-mer: " << maxKmer << " (count: " << maxCount << ")\n";
    
//     std::cout << "\nSample k-mers:\n";
//     int shown = 0;
//     for (const auto& [kmer, count] : results) {
//         if (shown++ >= 5) break;
//         std::cout << "  " << kmer << " -> " << count << "\n";
//     }

//     std::cout << "\nProcessing complete!\n";
//     return 0;
// }

// int main() {
//     const int k = 6;
//     const int m = 5;
//     const unsigned NUM_THREADS = 8;
//     const size_t HASH_TABLE_SIZE = 10'000'000;
//     const size_t MAX_PROBE_STEPS = 100;

//     using clock = std::chrono::high_resolution_clock;
//     auto total_start = clock::now();

//     // Step 1: Generate test FASTA
//     std::cout << "Generating test FASTA...\n";
//     generateTestFasta("test_large.fasta", 5'000'000);

//     // Step 2: Read bundles
//     std::cout << "Reading FASTA bundles...\n";
//     FastReader reader("GCA_000002315.5_GRCg6a_genomic.fna");
//     auto bundles = reader.readFile();
//     std::cout << "Read " << bundles.size() << " bundles\n";

//     // Step 3: Shared queue for KmerBlocks
//     std::queue<KmerBlock*> inputQueue;
//     std::mutex queueLock;
//     std::condition_variable cv;

//     // Step 4: Initialize Hasher with shared queue
//     std::cout << "Initializing Hasher...\n";
//     Hasher hasher(inputQueue, NUM_THREADS, HASH_TABLE_SIZE, MAX_PROBE_STEPS);

//     // Step 5: Launch worker threads to process queue
//     std::vector<std::thread> threads;
//     std::cout << "Launching " << NUM_THREADS << " worker threads...\n";
//     for (unsigned i = 0; i < NUM_THREADS; i++) {
//         threads.emplace_back(&Hasher::worker, &hasher, i);
//     }

//     // Step 6: Compute super-mers and push into queue (Phase 1 style)
//     std::cout << "Processing bundles and pushing super-mers to queue...\n";
//     std::atomic<size_t> bundleIndex(0);

//     // Parallelize super-mer computation per bundle
//     std::vector<std::thread> supermerThreads;
//     for (unsigned t = 0; t < NUM_THREADS; t++) {
//         supermerThreads.emplace_back([&]() {
//             while (true) {
//                 size_t idx = bundleIndex.fetch_add(1);
//                 if (idx >= bundles.size()) break;

//                 auto& b = bundles[idx];
//                 auto superMers = computeSuperMers(
//                     std::string(b.data.begin(), b.data.end()), m, k
//                 );

//                 for (auto& sm : superMers) {
//                     pushSuperMerToQueue(sm, k, inputQueue, queueLock, cv);
//                 }
//             }
//         });
//     }

//     // Wait for super-mer threads to finish
//     for (auto& t : supermerThreads) t.join();

//     std::cout << "All super-mers pushed to queue.\n";

//     // Step 7: Signal Hasher that no more data is coming
//     hasher.signalComplete();

//     // Step 8: Wait for Hasher threads to finish
//     for (auto& t : threads) t.join();

//     // Step 9: Merge results and output
//     std::cout << "Merging results from threads...\n";
//     hasher.mergeResults();
//     const auto& results = hasher.getResults();

//     std::cout << "Total unique k-mers: " << results.size() << "\n";
//     std::cout << "Writing results to output.txt...\n";
//     hasher.writeResults("output.txt");

//     // Step 10: Print stats
//     size_t maxCount = 0;
//     std::string maxKmer;
//     int shown = 0;
//     std::cout << "\n=== Statistics ===\n";
//     for (const auto& [kmer, count] : results) {
//         if (count > maxCount) {
//             maxCount = count;
//             maxKmer = kmer;
//         }
//         if (shown++ < 5) {
//             std::cout << "  " << kmer << " -> " << count << "\n";
//         }
//     }
//     std::cout << "Most frequent k-mer: " << maxKmer << " (count: " << maxCount << ")\n";

//     auto total_end = clock::now();
//     std::cout << "TOTAL runtime: "
//               << std::chrono::duration<double>(total_end - total_start).count()
//               << " sec\n";

//     std::cout << "Processing complete!\n";
//     return 0;
// }

int main() {
    const int k = 6;
    const int m = 5;
    const unsigned NUM_THREADS = 8;
    const size_t HASH_TABLE_SIZE = 10'000'000;
    const size_t MAX_PROBE_STEPS = 100;
    const size_t MAX_QUEUE_BATCH = 500; // batch size to push to queue

    using clock = std::chrono::high_resolution_clock;
    auto total_start = clock::now();

    std::cout << "Generating test FASTA...\n";
    generateTestFasta("test_large.fasta", 5'000'000);

    std::cout << "Reading FASTA bundles...\n";
    FastReader reader("GCA_000002315.5_GRCg6a_genomic.fna");
    auto bundles = reader.readFile();
    std::cout << "Read " << bundles.size() << " bundles\n";

    std::queue<KmerBlock*> inputQueue;
    std::mutex queueLock;
    std::condition_variable cv;

    std::cout << "Initializing Hasher...\n";
   Hasher hasher(inputQueue, NUM_THREADS, HASH_TABLE_SIZE, MAX_PROBE_STEPS);

    // Launch Hasher worker threads
    std::vector<std::thread> hashThreads;
    for (unsigned i = 0; i < NUM_THREADS; i++) {
        hashThreads.emplace_back(&Hasher::worker, &hasher, i);
    }

    std::cout << "Processing bundles and pushing super-mers to queue...\n";
    std::atomic<size_t> bundleIndex(0);

    // Producer threads: compute super-mers and push in batches
    std::vector<std::thread> producerThreads;
    for (unsigned t = 0; t < NUM_THREADS; t++) {
        producerThreads.emplace_back([&]() {
            std::vector<KmerBlock*> localBatch;
            while (true) {
                size_t idx = bundleIndex.fetch_add(1);
                if (idx >= bundles.size()) break;

                auto& b = bundles[idx];
                auto superMers = computeSuperMers(std::string(b.data.begin(), b.data.end()), m, k);

                for (auto& sm : superMers) {
                    size_t numKmers = sm.size() >= (size_t)k ? sm.size() - k + 1 : 0;
                    if (numKmers == 0) continue;

                    KmerBlock* block = new KmerBlock(numKmers);
                    for (size_t i = 0; i < numKmers; ++i)
                        block->kmers.push_back(sm.substr(i, k));

                    localBatch.push_back(block);

                    if (localBatch.size() >= MAX_QUEUE_BATCH) {
                        std::lock_guard<std::mutex> lock(queueLock);
                        for (auto blk : localBatch) inputQueue.push(blk);
                        cv.notify_all();
                        localBatch.clear();
                    }
                }
            }

            // Push any remaining batch
            if (!localBatch.empty()) {
                std::lock_guard<std::mutex> lock(queueLock);
                for (auto blk : localBatch) inputQueue.push(blk);
                cv.notify_all();
                localBatch.clear();
            }
        });
    }

    // Wait for producer threads to finish
    for (auto& t : producerThreads) t.join();

    std::cout << "All super-mers pushed to queue.\n";

    // Signal Hasher threads that no more data is coming
    hasher.signalComplete();

    // Wait for Hasher threads to finish
    for (auto& t : hashThreads) t.join();

    std::cout << "Merging results from threads...\n";
    hasher.mergeResults();
    const auto& results = hasher.getResults();

    std::cout << "Total unique k-mers: " << results.size() << "\n";
    hasher.writeResults("output.txt");

    size_t maxCount = 0;
    std::string maxKmer;
    int shown = 0;
    std::cout << "\n=== Statistics ===\n";
    for (const auto& [kmer, count] : results) {
        if (count > maxCount) { maxCount = count; maxKmer = kmer; }
        if (shown++ < 5) std::cout << "  " << kmer << " -> " << count << "\n";
    }
    std::cout << "Most frequent k-mer: " << maxKmer << " (count: " << maxCount << ")\n";

    auto total_end = clock::now();
    std::cout << "TOTAL runtime: "
              << std::chrono::duration<double>(total_end - total_start).count()
              << " sec\n";

    std::cout << "Processing complete!\n";
    return 0;
}
