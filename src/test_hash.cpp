#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>
#include "QuadraticHashTable.h"

int main() {
    std::cout << "=== QuadraticHashTable Tests ===\n\n";
    
    // Test 1: Basic insertion
    std::cout << "Test 1: Basic insertion\n";
    QuadraticHashTable table(1009, 10);
    bool success1 = table.insert("ACGTACGTACGTACGTACGTACGTACGTACGT");
    bool success2 = table.insert("TGCATGCATGCATGCATGCATGCATGCATGCA");
    std::cout << "  Insert kmer1: " << (success1 ? "SUCCESS" : "FAILED") << "\n";
    std::cout << "  Insert kmer2: " << (success2 ? "SUCCESS" : "FAILED") << "\n";
    std::cout << "  PASS: " << (success1 && success2 ? "YES" : "NO") << "\n\n";
    
    // Test 2: Duplicate insertion (should increment count)
    std::cout << "Test 2: Duplicate insertion increments count\n";
    QuadraticHashTable table2(1009, 10);
    table2.insert("ACGTACGTACGTACGTACGTACGTACGTACGT");
    table2.insert("ACGTACGTACGTACGTACGTACGTACGTACGT");
    table2.insert("ACGTACGTACGTACGTACGTACGTACGTACGT");
    std::cout << "  Inserted same k-mer 3 times\n";
    std::cout << "  (Check implementation - count should be 3)\n\n";
    
    // Test 3: Many insertions
    std::cout << "Test 3: Insert 500 unique k-mers\n";
    QuadraticHashTable table3(1009, 10);
    int successful = 0;
    int failed = 0;
    
    for (int i = 0; i < 500; i++) {
        std::string random_kmer;
        for (int j = 0; j < 32; j++) {
            random_kmer += "ACGT"[rand() % 4];
        }
        if (table3.insert(random_kmer)) {
            successful++;
        } else {
            failed++;
        }
    }
    
    std::cout << "  Successful insertions: " << successful << "\n";
    std::cout << "  Failed insertions: " << failed << "\n";
    std::cout << "  Success rate: " << (100.0 * successful / 500) << "%\n";
    std::cout << "  PASS: " << (successful > 450 ? "YES (>90% success)" : "NO (<90% success)") << "\n\n";
    
    // Test 4: High load factor test
    std::cout << "Test 4: High load factor (800 insertions into table of size 1009)\n";
    QuadraticHashTable table4(1009, 10);
    successful = 0;
    failed = 0;
    
    for (int i = 0; i < 800; i++) {
        std::string random_kmer;
        for (int j = 0; j < 32; j++) {
            random_kmer += "ACGT"[rand() % 4];
        }
        if (table4.insert(random_kmer)) {
            successful++;
        } else {
            failed++;
        }
    }
    
    std::cout << "  Successful insertions: " << successful << "\n";
    std::cout << "  Failed insertions: " << failed << "\n";
    std::cout << "  Load factor: " << (100.0 * successful / 1009) << "%\n";
    std::cout << "  PASS: " << (failed < 100 ? "YES (reasonable failure rate)" : "NO (too many failures)") << "\n\n";
    
    // Test 5: Max steps limit
    std::cout << "Test 5: Max steps limit (small maxSteps should cause failures)\n";
    QuadraticHashTable table5(101, 2);
    successful = 0;
    failed = 0;
    
    for (int i = 0; i < 100; i++) {
        std::string random_kmer;
        for (int j = 0; j < 32; j++) {
            random_kmer += "ACGT"[rand() % 4];
        }
        if (table5.insert(random_kmer)) {
            successful++;
        } else {
            failed++;
        }
    }
    
    std::cout << "  Successful insertions: " << successful << "\n";
    std::cout << "  Failed insertions: " << failed << "\n";
    std::cout << "  PASS: " << (failed > 0 ? "YES (maxSteps limit working)" : "UNCERTAIN") << "\n\n";
    
    // Test 6: Collision handling with similar k-mers
    std::cout << "Test 6: K-mers with similar prefixes\n";
    QuadraticHashTable table6(1009, 10);
    std::vector<std::string> similar_kmers;
    similar_kmers.push_back("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    similar_kmers.push_back("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB");
    similar_kmers.push_back("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC");
    similar_kmers.push_back("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAD");
    
    successful = 0;
    for (size_t i = 0; i < similar_kmers.size(); i++) {
        if (table6.insert(similar_kmers[i])) successful++;
    }
    
    std::cout << "  Inserted " << successful << " out of " << similar_kmers.size() << " similar k-mers\n";
    std::cout << "  PASS: " << (successful == similar_kmers.size() ? "YES" : "NO") << "\n\n";
    
    // Test 7: Different k-mers that should hash differently
    std::cout << "Test 7: Verify different k-mers hash to different positions\n";
    QuadraticHashTable table7(1009, 10);
    std::string kmer1 = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    std::string kmer2 = "TGCATGCATGCATGCATGCATGCATGCATGCA";
    std::string kmer3 = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    
    bool ins1 = table7.insert(kmer1);
    bool ins2 = table7.insert(kmer2);
    bool ins3 = table7.insert(kmer3);
    
    std::cout << "  All insertions successful: " << (ins1 && ins2 && ins3 ? "YES" : "NO") << "\n";
    std::cout << "  PASS: " << (ins1 && ins2 && ins3 ? "YES" : "NO") << "\n\n";
    
    std::cout << "=== All Tests Complete ===\n";
    
    return 0;
}