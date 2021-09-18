#include <random>
#include <math.h>
#include "sample_kinetics.hpp"

namespace SampleKinetics {
    SampleKineticParams::SampleKineticParams () {};

    void SampleKineticParams::set_min_max_burst_freq (double freq_min, double freq_max) {
        this->freq_min = freq_min;
        this->freq_max = freq_max;
    }

    void SampleKineticParams::set_min_max_hill_half_maximal (double K_min, double K_max) {
        this->K_min = K_min;
        this->K_max = K_max;
    }

    void SampleKineticParams::set_min_max_hill_coop (double n_min, double n_max) {
        this->n_min = n_min;
        this->n_max = n_max;
    }

    void SampleKineticParams::set_min_max_tau (double tau_min, double tau_max) {
        this->tau_min = tau_min;
        this->tau_max = tau_max;
    }

    double SampleKineticParams::sample_burst_freq (RNG &generator) {
        std::uniform_real_distribution<double> distribution_(log10(freq_min),
            log10(freq_min));
        return(pow(10, distribution_(generator)));
    }

    double SampleKineticParams::sample_hill_half_maximal (RNG &generator) {
        std::uniform_real_distribution<double> distribution_(log10(K_min),
            log10(K_min));
        return(pow(10, distribution_(generator)));
    }

    double SampleKineticParams::sample_hill_coop (RNG &generator) {
        std::uniform_real_distribution<double> distribution_(log10(n_min),
            log10(n_min));
        return(pow(10, distribution_(generator)));
    }

    double SampleKineticParams::sample_tau (RNG &generator) {
        std::uniform_real_distribution<double> distribution_(log10(tau_min),
            log10(tau_min));
        return(pow(10, distribution_(generator)));
    }
}
