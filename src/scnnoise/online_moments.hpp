/* Copyright (C) 2021 Tarun Mahajan - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GNU AGPLv3 license.
 *
 * You should have received a copy of the GNU AGPLv3 license with
 * this file. If not, please write to: tarunm3@illinois.edu, or 
 * visit : https://www.gnu.org/licenses/agpl-3.0.en.html
 */

// Class for computing online mean and variance
#ifndef ONLINEMOMENTS_H
#define ONLINEMOMENTS_H

#include <vector>

namespace OnlineMomentsSpace {
    struct aggregate_struct {
        double mean_;
        double M2;
        double count;
    };

    class OnlineMoments {
    private:

    public:
        /*Data member */
        aggregate_struct existing_aggregate;
        /* Member functions */
        // Constructor
        OnlineMoments (std::vector<double> samples);

        void update_moments (double new_sample);
    };
}

#endif
