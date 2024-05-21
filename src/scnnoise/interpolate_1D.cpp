/* Copyright (C) 2021 Tarun Mahajan - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GNU AGPLv3 license.
 *
 * You should have received a copy of the GNU AGPLv3 license with
 * this file. If not, please write to: tarunm3@illinois.edu, or 
 * visit : https://www.gnu.org/licenses/agpl-3.0.en.html
 */

// Class for computing online mean and variance
#include <vector>
#include <algorithm>
#include "interpolate_1D.hpp"
#include <math.h>

namespace LinearInterpolate {
    LinearInterpolation_1D::LinearInterpolation_1D (std::vector<double> x_samples, std::vector<double> y_samples) {
        data_pairs.x_vals = x_samples;
        data_pairs.y_vals = y_samples;
    }

    // Interpolate single sample
    double LinearInterpolation_1D::interpolate_1D (double x_sample, double x_lower, double x_upper,
        double y_lower, double y_upper) {
            double y_sample = (y_upper - y_lower);
            y_sample /= (x_upper - x_lower);
            y_sample *= (x_sample - x_lower);
            y_sample += y_lower;
        }

    // Interpolate values
    std::vector<double> LinearInterpolation_1D::interpolate_samples (std::vector<double> &x_samples) {
        std::vector<double> y_samples(x_samples.size(), -1);
        // std::vector<int>::iterator low,up;
        for (std::size_t b = 0; b < x_samples.size(); ++b) {
            auto up = std::upper_bound (data_pairs.x_vals.begin(),
                data_pairs.x_vals.end(), x_samples[b]);
            if (up != data_pairs.x_vals.end() && up - data_pairs.x_vals.begin() > 0) {
                unsigned int id_lower = up - data_pairs.x_vals.begin() - 1;
                if (id_lower > 0) {
                    if (data_pairs.x_vals[id_lower] == x_samples[b]) {
                        id_lower -= 1;
                    }
                    y_samples[b] = interpolate_1D (x_samples[b],
                        data_pairs.x_vals[id_lower],
                        data_pairs.x_vals[up - data_pairs.x_vals.begin()],
                        data_pairs.y_vals[id_lower],
                        data_pairs.y_vals[up - data_pairs.x_vals.begin()]);
                }

            }
        }
        return y_samples;
    }
}
