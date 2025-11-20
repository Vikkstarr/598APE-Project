#include <string>
#include <vector>

std::vector<std::string> generateKmers(const std::string& seq, int k);
std::string computeMinimizer(const std::string& kmer, int m, int k);
std::vector<std::string> computeAllMinimizers(const std::string& seq, int m, int k);
std::vector<std::string> computeSuperMers(const std::string& seq, int m, int k);

// add the method header for superMerToKmers
std::vector<std::string> superMerToKmers(const std::string& superMer, int k);

