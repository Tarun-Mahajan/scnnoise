#include <random>
#include <math.h>
#include "sample_kinetics.hpp"

namespace SampleKinetics {
    SampleKineticParams::SampleKineticParams () {};

    void SampleKineticParams::set_min_max_burst_freq (double freq_min, double freq_max) {
        this->freq_min = freq_min;
        this->freq_max = freq_max;
    }

    void SampleKineticParams::set_min_max_burst_size (int burst_size_min, int burst_size_max) {
        this->burst_size_min = burst_size_min;
        this->burst_size_max = burst_size_max;
    }

    void SampleKineticParams::set_min_max_hill_half_maximal (double K_min, double K_max) {
        this->K_min = K_min;
        this->K_max = K_max;
    }

    void SampleKineticParams::set_min_max_hill_coop (int n_min, int n_max) {
        this->n_min = n_min;
        this->n_max = n_max;
    }

    void SampleKineticParams::set_min_max_tau (double tau_min, double tau_max) {
        this->tau_min = tau_min;
        this->tau_max = tau_max;
    }

    double SampleKineticParams::sample_burst_freq (RNG &generator) {
        std::uniform_real_distribution<double> distribution_(log10(freq_min),
            log10(freq_max));
        return(pow(10, distribution_(generator)));
    }

    double SampleKineticParams::sample_burst_size (RNG &generator) {
        std::uniform_int_distribution<int> distribution_(burst_size_min,
            burst_size_max);
        int burst_ = distribution_(generator);
        return(double(burst_));
    }

    double SampleKineticParams::sample_hill_half_maximal (RNG &generator) {
        std::uniform_real_distribution<double> distribution_(log10(K_min),
            log10(K_max));
        return(pow(10, distribution_(generator)));
    }

    double SampleKineticParams::sample_hill_coop (RNG &generator) {
        // std::uniform_real_distribution<double> distribution_(log10(n_min),
        //     log10(n_max));
        // return(pow(10, distribution_(generator)));

        std::uniform_int_distribution<int> distribution_(n_min,
            n_max);
        return(double(distribution_(generator)))
    }

    double SampleKineticParams::sample_tau (RNG &generator) {
        std::uniform_real_distribution<double> distribution_(log10(tau_min),
            log10(tau_max));
        return(pow(10, distribution_(generator)));
    }

    bool SampleKineticParams::sample_activator_status (RNG &generator) {
        std::uniform_int_distribution<int> distribution_(0, 1);
        if (distribution_(generator) == 0) {
            return(false);
        }else{
            return(true);
        }
    }
}
