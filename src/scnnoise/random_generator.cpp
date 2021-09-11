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
