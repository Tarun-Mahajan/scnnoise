// Class for computing online mean and variance
#include <vector>
#include "interpolate_1D.hpp"
#include <math.h>

namespace LinearInterpolate {
    LinearInterpolation_1D::LinearInterpolation_1D (std::vector<double> x_samples, std::vector<double> y_samples) {
        data_struct.x_vals = x_samples;
        data_struct.y_vals = y_samples;
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
        std::vector<int>::iterator low,up;
        for (std::size_t b = 0; b < x_samples.size(); ++b) {
            up = std::upper_bound (data_struct.x_vals.begin(),
                data_struct.x_vals.end(), x_samples[b]);
            if (up != data_struct.x_vals.end() && up - data_struct.x_vals.begin() > 0) {
                unsigned int id_lower = up - data_struct.x_vals.begin() - 1;
                if (id_lower > 0) {
                    if (data_struct.x_vals[id_lower] == x_samples[b]) {
                        id_lower -= 1;
                    }
                    y_samples[b] = interpolate_1D (x_samples[b],
                        data_struct.x_vals[id_lower],
                        data_struct.x_vals[up - data_struct.x_vals.begin()],
                        data_struct.y_vals[id_lower],
                        data_struct.y_vals[up - data_struct.x_vals.begin()]);
                }

            }
        }
        return y_samples;
    }
}
