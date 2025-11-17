#ifndef QUADRATIC_HASH_TABLE_H
#define QUADRATIC_HASH_TABLE_H

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>

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

                for (int i = 0; i < tableSize; i++) {
                    keys[i] = "";
                    values[i] = 0;
                }
        }

        bool insert(const std::string& kmer) {
            size_t i = 0;
            size_t hashPos;
            
            while (true) {
                hashPos = hash<32>(kmer, i) % tableSize;

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
            const uint64_t* data = reinterpret_cast<const uint64_t*>(x.data());
            uint64_t res = 0;

            for(int c = 0; c < (k + 7) / 8; ++c) {
                res = 453569 * res + data[c];
            }

            return res + 5696063 * i * i;
        }

        void exportToMap(std::unordered_map<std::string, size_t>& map) const {
            for (size_t i = 0; i < tableSize; i++) {
                if (keys[i] != "") {
                    map[keys[i]] += values[i];
                }
            }
        }
};

#endif