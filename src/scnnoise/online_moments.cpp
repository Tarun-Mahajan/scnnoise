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
#include "online_moments.hpp"
#include <math.h>

namespace OnlineMomentsSpace {
    OnlineMoments::OnlineMoments (std::vector<double> samples) {
        existing_aggregate.mean_ = 0;
        existing_aggregate.M2 = 0;
        existing_aggregate.count = samples.size();

        // compute mean
        for (std::size_t i = 0; i < samples.size(); ++i) {
            existing_aggregate.mean_ += samples[i];
        }
        existing_aggregate.mean_ /= samples.size();

        // compute sum squared difference from mean
        for (std::size_t i = 0; i < samples.size(); ++i) {
            existing_aggregate.M2 += pow(samples[i] - existing_aggregate.mean_, 2.0);
        }
    }

    // Implement Welford's algorithm
    void OnlineMoments::update_moments (double new_sample) {
        existing_aggregate.count += 1;
        double delta = new_sample - existing_aggregate.mean_;
        existing_aggregate.mean_ += delta / existing_aggregate.count;
        double delta2 = new_sample - existing_aggregate.mean_;
        existing_aggregate.M2 += delta * delta2;
    }
}
