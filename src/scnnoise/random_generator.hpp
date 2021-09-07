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
