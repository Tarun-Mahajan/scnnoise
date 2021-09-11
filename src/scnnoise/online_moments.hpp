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
