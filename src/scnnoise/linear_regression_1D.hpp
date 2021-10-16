// Graph header file
#ifndef LINEAR_REGRESS_H
#define LINEAR_REGRESS_H

#include <vector>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

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
