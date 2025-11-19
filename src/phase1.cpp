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
    //#pragma omp parallel for
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

// lowkey same as readFile but pushing to queue in batches of 10
// with locking and cv notifying

void pushSuperMersToQueue(const std::vector<std::string>& superMers, int k,
                           std::queue<KmerBlock*>& inputQueue,
                           std::mutex& queueLock,
                           std::condition_variable& cv) {
    // first u want to reserve a batching of 10
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
            // push to shared queue
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
    // push any remaining blocks
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

 #include <fstream>
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

    // 1. Generate ASTA file
    std::cout << "Generating test FASTA...\n";
    generateTestFasta("test_large.fasta", 5'000'000);

    // 2. Read the FASTA using your FastReader
    std::cout << "Reading FASTA bundles...\n";
    FastReader reader("test_large.fasta");
    auto bundles = reader.readFile(); // 1 MB bundles
    std::cout << "Read " << bundles.size() << " bundles\n";

    // 3. Compute super-mers for each bundle
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

    // 4. Build input queue with kmer blocks
    std::queue<KmerBlock*> inputQueue;
    std::mutex queueLock;
    std::condition_variable cv;

    // 5. Push super-mers into queue in batches
    std::cout << "Pushing super-mers to queue...\n";
    pushSuperMersToQueue(superMers, k, inputQueue, queueLock, cv);

    std::cout << "Queue size after batching: " << inputQueue.size() << "\n";

    // 6. Inspect a few blocks
    int inspectN = 3;
    for (int i = 0; i < inspectN && !inputQueue.empty(); i++) {
        KmerBlock* block = inputQueue.front();
        inputQueue.pop();

        std::cout << "Block " << i
                  << " has " << block->kmers.size()
                  << " kmers\n";

        // Print the first few kmers from this block
        for (int j = 0; j < 5 && j < (int)block->kmers.size(); j++) {
            std::cout << "  kmer[" << j << "] = " << block->kmers[j] << "\n";
        }
    }

    std::cout << "Queue testing complete.\n";
    return 0;
}
