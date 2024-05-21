/* Copyright (C) 2021 Tarun Mahajan - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GNU AGPLv3 license.
 *
 * You should have received a copy of the GNU AGPLv3 license with
 * this file. If not, please write to: tarunm3@illinois.edu, or 
 * visit : https://www.gnu.org/licenses/agpl-3.0.en.html
 */

// scNNoiSE header file
#ifndef SAMPLE_KINETICS_H
#define SAMPLE_KINETICS_H

namespace SampleKinetics {
    class SampleKineticParams {
    private:
        /* data */
        double freq_min;
        double freq_max;
        double K_min;
        double K_max;
        double n_min;
        double n_max;
        double tau_min;
        double tau_max;
        int burst_size_min;
        int burst_size_max;


    public:
        typedef std::mt19937 RNG;
        /* Member functions */
        // Constructor
        SampleKineticParams ();

        // Set min, max for burst freq.
        void set_min_max_burst_freq (double freq_min, double freq_max);

        // Set min, max for burst size
        void set_min_max_burst_size (int burst_size_min, int burst_size_max);

        // Set min, max for hill half maximal
        void set_min_max_hill_half_maximal (double K_min, double K_max);

        // Set min, max for hill cooperativity
        void set_min_max_hill_coop (int n_min, int n_max);

        // Set min, max for tau (1/degradation_rate)
        void set_min_max_tau (double tau_min, double tau_max);

        // Sample burst freq
        double sample_burst_freq (RNG &generator);

        // Sample burst size
        double sample_burst_size (RNG &generator);

        // Sample hill half maximal
        double sample_hill_half_maximal (RNG &generator);

        // Sample hill cooperativity
        double sample_hill_coop (RNG &generator);

        // Sample tau (1/degradation_rate)
        double sample_tau (RNG &generator);

        // sample activator status
        bool sample_activator_status (RNG &generator);


    };
}
#endif
