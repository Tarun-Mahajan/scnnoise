/* Copyright (C) 2021 Tarun Mahajan - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GNU AGPLv3 license.
 *
 * You should have received a copy of the GNU AGPLv3 license with
 * this file. If not, please write to: tarunm3@illinois.edu, or 
 * visit : https://www.gnu.org/licenses/agpl-3.0.en.html
 */

// Graph header file
#ifndef LINEAR_REGRESS_H
#define LINEAR_REGRESS_H

#include <vector>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>

namespace LinearRegressionSpace {
    class LinearRegression1D {
    public:
        // data members
        unsigned int num_samples;
        std::vector<double> input_samples;
        std::vector<double> output_samples;
        double slope;
        double intercept;
        double t_statistics_;
        double alpha_;
        double p_value;
        /* Memeber functions */
        // Constructor
        LinearRegression1D ();

        void set_sample_size (unsigned int num_samples);

        void set_input_samples (std::vector<double> input_samples);

        void set_output_samples (std::vector<double> output_samples);

        void set_alpha (double alpha_);

        void estimation_ ();

        void comput_p_value ();

        void t_test_ ();

        bool reached_steady_state ();
    };
}

#endif
