/* Copyright (C) 2021 Tarun Mahajan - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GNU AGPLv3 license.
 *
 * You should have received a copy of the GNU AGPLv3 license with
 * this file. If not, please write to: tarunm3@illinois.edu, or 
 * visit : https://www.gnu.org/licenses/agpl-3.0.en.html
 */


// Solve for mean expression using hill function regulation
#ifndef FIND_MEAN_H
#define FIND_MEAN_H

#include <vector>

namespace ScnnoiseInterface {
    class CostFunctionMean {
    public:
        virtual bool Evaluate(double const* const* parameters,
                            double* residuals,
                            double** jacobians) = 0;
        const vector<int32>& parameter_block_sizes();
        int num_residuals() const;

    protected:
        vector<int32>* mutable_parameter_block_sizes();
        void set_num_residuals(int num_residuals);
    };

    class LinearInterpolation_1D {
    private:

    public:
        /*Data member */
        data_struct data_pairs;
        /* Member functions */
        // Constructor
        LinearInterpolation_1D (std::vector<double> x_samples,
            std::vector<double> y_samples);

        double interpolate_1D (double x_sample, double x_lower, double x_upper,
            double y_lower, double y_upper);

        std::vector<double> interpolate_samples (std::vector<double> &x_samples);
    };
}

#endif
