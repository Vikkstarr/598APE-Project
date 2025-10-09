#include <iostream>
#include <fstream>
#include <vector>
#include <string>

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
        std::ifstream in(path, std::ios::binary);
        if (!in) {
            throw std::runtime_error("Could not open file: " + path);
        }

        std::vector<FastBundle> bundles;
        std::vector<char> buffer(blockSize);

        while (in) {
            in.read(buffer.data(), blockSize);
            std::streamsize n = in.gcount();
            if (n > 0) {
                FastBundle bundle(blockSize);
                bundle.addBlock(buffer.data(), n);
                bundle.finalize();
                bundles.push_back(std::move(bundle));
            }
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

// super-mer of a genome sequence is defined as a substring of maximal 
// length such that all k-mers on that substring share the same minimizer. 
// Hereby, a minimizer of a k-mer is defined as a substring of fixed length 
// m that is minimal with respect to some total ordering on strings of length 
// m. Thus, contiguous k-mers of a genome read are joined to a super-mer if they 
// share the same minimizer.”


std::string computeMinimizer(const std::string &seq, int m, int k) {
    // so basically take the first m-mer of the first k-mer and itetate thorugh every one

    if ((int)seq.size() < k || m > k) return "";  // too short or invalid

    std::string minimizer = seq.substr(0, m);
    for (size_t i = 1; i <= seq.size() - m; i++) {
        std::string current_mmer = seq.substr(i, m);
        if (current_mmer < minimizer) {
            minimizer = current_mmer;
        }
    }
}

std::vector<std::string> computeAllMinimizers(const std::string &seq, int m, int k) {
    if ((int)seq.size() < k || m > k) return;  // too short or invalid
    auto kmers = generateKmers(seq, k);

    std::vector<std::string> minimizers;
    for (const auto &kmer : kmers) {
        minimizers.push_back(computeMinimizer(kmer, m, k));
    }
    return minimizers;
}

// i have to now group super mers based on minimizers


int main() {
    FastReader reader("example_sequences.fasta");  // pick one FASTQ/FASTA file
    auto bundles = reader.readFile();

    std::cout << "Read " << bundles.size() << " bundles\n";
    for (size_t i = 0; i < bundles.size(); i++) {
        std::cout << "Bundle " << i << " size = " << bundles[i].data.size() << " bytes\n";
    }
}
