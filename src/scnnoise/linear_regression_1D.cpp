#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#define _USE_MATH_DEFINES
#include <math.h>
#include <boost/math/distributions/students_t.hpp>
#include "linear_regression_1D.hpp"

using namespace boost::math;

namespace LinearRegressionSpace {
    LinearRegression1D::LinearRegression1D () {
        alpha_ = 0.001;
    };

    void LinearRegression1D::set_sample_size (unsigned int num_samples) {
        this->num_samples = num_samples;
    }

    void LinearRegression1D::set_input_samples (std::vector<double> input_samples) {
        this->input_samples = input_samples;
    }

    void LinearRegression1D::set_output_samples (std::vector<double> output_samples) {
        this->output_samples = output_samples;
    }

    void LinearRegression1D::set_alpha (double alpha_) {
        this->alpha_ = alpha_;
    }

    void LinearRegression1D::estimation_ () {

        // XTX data covariance matrix
        std::vector<std::vector<double>> XTX(2);
        std::vector<std::vector<double>> XTX_inv(2);
        for (unsigned int i = 0; i < 2; ++i) {
            XTX[i].resize(2);
            XTX_inv[i].resize(2);
        }
        XTX[0][0] = num_samples;
        XTX[0][1] = 0;
        XTX[1][0] = XTX[0][1];
        XTX[1][1] = 0;
        for (std::size_t i = 0; i < output_samples.size(); ++i) {
            XTX[1][1] += pow(input_samples[i], 2);
            XTX[0][1] += input_samples[i];
        }
        XTX[1][0] = XTX[0][1];

        // determinant for XTX data covariance matrix
        double delta_XTX = XTX[0][0] * XTX[1][1] - XTX[0][1] * XTX[1][0];

        // Inverse for XTX data covariance matrix
        XTX_inv[0][0] = XTX[1][1]/delta_XTX;
        // std::cout << "Error here = " << delta_XTX << std::endl;
        XTX_inv[0][1] = -XTX[1][0]/delta_XTX;
        XTX_inv[1][0] = -XTX[0][1]/delta_XTX;
        XTX_inv[1][1] = XTX[0][0]/delta_XTX;


        // XTy data output dot product
        std::vector<double> XTy(2, 0);
        XTy[0] = std::accumulate(output_samples.begin(),
            output_samples.end(),
            decltype(output_samples)::value_type(0));
        XTy[1] = 0;
        for (std::size_t i = 0; i < output_samples.size(); ++i) {
            XTy[1] += input_samples[i] * output_samples[i];
        }

        // Estimate intercept and slope
        intercept = XTX_inv[0][0] * XTy[0] + XTX_inv[0][1] * XTy[1];
        slope = XTX_inv[1][0] * XTy[0] + XTX_inv[1][1] * XTy[1];
    }

    void LinearRegression1D::comput_p_value () {
        students_t dist_(num_samples - 2);
        p_value = cdf(complement(dist_, fabs(t_statistics_)));
    }

    void LinearRegression1D::t_test_ () {
        std::vector<double> output_predicted(num_samples, 0);
        double input_mean = 0;
        double SE_slope = 0;
        double SE_num = 0;
        double SE_denom = 0;

        for (std::size_t i = 0; i < num_samples; ++i) {
            input_mean += input_samples[i];
        }
        input_mean /= num_samples;

        for (std::size_t i = 0; i < num_samples; ++i) {
            output_predicted[i] = input_samples[i] * slope + intercept;
        }

        for (std::size_t i = 0; i < num_samples; ++i) {
            SE_num += pow(output_samples[i] - output_predicted[i], 2);
            SE_denom += pow(input_samples[i] - input_mean, 2);
        }
        SE_num /= num_samples - 2;
        SE_num = sqrt(SE_num);
        SE_denom = sqrt(SE_denom);

        SE_slope = SE_num/SE_denom;
        t_statistics_ = slope/SE_slope;
        comput_p_value();
    }

    bool LinearRegression1D::reached_steady_state () {
        estimation_();
        t_test_();
        if (p_value > alpha_/2) {
            return false;
        }else{
            return true;
        }
    }
}
