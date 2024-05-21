/* Copyright (C) 2021 Tarun Mahajan - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GNU AGPLv3 license.
 *
 * You should have received a copy of the GNU AGPLv3 license with
 * this file. If not, please write to: tarunm3@illinois.edu, or 
 * visit : https://www.gnu.org/licenses/agpl-3.0.en.html
 */

#include "random_generator.hpp"

namespace RandomGenerator {
    RandomGeneratorSeeds::RandomGeneratorSeeds (unsigned int num_seeds) {
        this->num_seeds = num_seeds;
    }

    std::vector<std::uint_least32_t> RandomGeneratorSeeds::generate_ (std::seed_seq &sd) {
        std::vector<std::uint_least32_t> rd_seeds(num_seeds);
        sd.generate(rd_seeds.begin(), rd_seeds.end());
        return rd_seeds;
    }
}
