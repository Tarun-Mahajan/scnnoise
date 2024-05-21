/* Copyright (C) 2021 Tarun Mahajan - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GNU AGPLv3 license.
 *
 * You should have received a copy of the GNU AGPLv3 license with
 * this file. If not, please write to: tarunm3@illinois.edu, or 
 * visit : https://www.gnu.org/licenses/agpl-3.0.en.html
 */

// Generate seeds for random number generator
#ifndef RANDOM_GEN
#define RANDOM_GEN

#include <random>

namespace RandomGenerator {
    class RandomGeneratorSeeds {
    private:

    public:
        /*Data member */
        unsigned int num_seeds;
        /* Member functions */
        // Constructor
        RandomGeneratorSeeds (unsigned int num_seeds);

        std::vector<std::uint_least32_t> generate_ (std::seed_seq &sd);
    };
}

#endif
