#include <iostream>
#include <fstream>
#include "QuadraticHashTable.h"

// Simple test to verify duplicate k-mer counting works
int main() {
    QuadraticHashTable table(1000, 100);
    
    std::cout << "Testing duplicate k-mer counting...\n\n";
    
    // Insert the same k-mer 5 times
    std::string kmer1 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";  // 31 A's
    
    std::cout << "Inserting '" << kmer1 << "' 5 times...\n";
    for (int i = 0; i < 5; i++) {
        bool success = table.insert(kmer1);
        std::cout << "  Insert " << (i+1) << ": " << (success ? "success" : "failed") << "\n";
    }
    
    // Insert a different k-mer 3 times
    std::string kmer2 = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";  // 31 T's
    
    std::cout << "\nInserting '" << kmer2 << "' 3 times...\n";
    for (int i = 0; i < 3; i++) {
        bool success = table.insert(kmer2);
        std::cout << "  Insert " << (i+1) << ": " << (success ? "success" : "failed") << "\n";
    }
    
    // Export to map and check counts
    std::unordered_map<std::string, size_t> map;
    table.exportToMap(map);
    
    std::cout << "\nResults:\n";
    for (const auto& [kmer, count] : map) {
        std::cout << "  " << kmer << " -> " << count << "\n";
    }
    
    // Verify
    std::cout << "\n";
    if (map[kmer1] == 5) {
        std::cout << "✓ PASS: kmer1 counted correctly (5)\n";
    } else {
        std::cout << "✗ FAIL: kmer1 count = " << map[kmer1] << ", expected 5\n";
    }
    
    if (map[kmer2] == 3) {
        std::cout << "✓ PASS: kmer2 counted correctly (3)\n";
    } else {
        std::cout << "✗ FAIL: kmer2 count = " << map[kmer2] << ", expected 3\n";
    }
    
    return 0;
}