#ifndef QUADRATIC_HASH_TABLE_H
#define QUADRATIC_HASH_TABLE_H

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <functional>
#include <iostream>

class QuadraticHashTable {
    private:
        std::vector<std::string> keys;
        std::vector<size_t> values;
        size_t tableSize;
        size_t numElements;
        size_t maxSteps;

    public:
        QuadraticHashTable(size_t size = 1009, size_t maxSteps = 5) 
            : tableSize(size), numElements(0), maxSteps(maxSteps) {
                keys.resize(size);
                values.resize(size);

                for (size_t i = 0; i < tableSize; i++) {
                    keys[i] = "";
                    values[i] = 0;
                }
        }

        bool insert(const std::string& kmer) {
            size_t i = 0;
            size_t hashPos;
            
            while (true) {
                hashPos = computeHash(kmer, i) % tableSize;

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

        uint64_t computeHash(const std::string& str, size_t i) const {
            std::hash<std::string> hasher;
            uint64_t baseHash = hasher(str);
            
            return baseHash + 5696063 * i * i;
        }

        void exportToMap(std::unordered_map<std::string, size_t>& map) const {
            for (size_t i = 0; i < tableSize; i++) {
                if (keys[i] != "") {
                    map[keys[i]] += values[i];
                }
            }
        }
        
        void printStats() const {
            size_t occupied = 0;
            size_t totalCount = 0;
            for (size_t i = 0; i < tableSize; i++) {
                if (keys[i] != "") {
                    occupied++;
                    totalCount += values[i];
                }
            }
            std::cout << "  Slots occupied: " << occupied << "/" << tableSize << "\n";
            std::cout << "  Total k-mer count: " << totalCount << "\n";
        }
};

#endif