/* Copyright (C) 2021 Tarun Mahajan - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GNU AGPLv3 license.
 *
 * You should have received a copy of the GNU AGPLv3 license with
 * this file. If not, please write to: tarunm3@illinois.edu, or 
 * visit : https://www.gnu.org/licenses/agpl-3.0.en.html
 */

// Optimized gillespie header file
#ifndef SDMCellCycle_H
#define SDMCellCycle_H

#include <vector>
#include <string>
#include "scnnoise.hpp"
#include "gillespieSSA.hpp"

namespace ScnnoiseInterface {
    struct generate_cell_cycle_time_struct {
        std::string mechanism;
        std::string distribution_;
        std::map<std::string, double> distribution_params;
        double replication_time_factor;
    };

    const std::map<std::string, double> distribution_params_init =
        {{"alpha", 12.0}, {"beta", 1.0/0.0083}, {"multiplicative_factor", 1.0}};
    /********************************************//**
    \brief A class for Gillespie's stochastic simulation algorithm.

    The Gillespie class creates an object to perform
    exact stochastic simulation for any given chemical
    reaction network.
    ***********************************************/

    class gillespieSDMCellCycle : public GillespieSSA {
    private:

    public:
        /*Data member */
        std::string current_cell_cycle_state;
        double cell_cycle_start_time;
        double cell_cycle_length;
        double cur_time;
        double next_time;
        generate_cell_cycle_time_struct cell_cycle_params;
        std::vector<double> dosage_compensation;
        double replication_time_factor;
        bool is_frozen_cell_cycle;
        std::string frozen_state;
        /* Member functions */
        // Constructor
        gillespieSDMCellCycle (int num_genes, std::string gene_filepath,
            std::string molecule_count_filepath,
            std::string count_save_file, bool keep_GRN,
            std::string GRN_filepath, int num_timepoints_save);

        void sort_reaction (int &rxn_selected) override;
        void update_cell_cycle_state (double next_time, double cur_time,
            RNG& generator) override;

        // distribution_params_init["alpha"] = 12;
        // distribution_params_init["beta"] = 1.0/0.0083;
        double generate_cell_cycle_time ();
        void set_cell_cycle_params (std::string mechanism = "timing",
            std::string distribution_ = "gamma",
            std::map<std::string, double> distribution_params = distribution_params_init,
            double replication_time_factor = 0.5);

        void set_dosage_compensation (
            std::vector<double> dosage_compensation);

        void init_cell_cycle_state (RNG &generator, double cur_time, 
            bool init_dosage_comp_adj=true) override;
        void cell_division (RNG& generator);
        void remove_dosage_compensation ();
        void update_propensity_cell_cycle ();
        void update_propensity_cell_division ();
        void sample_cell_cycle_time (RNG &generator);
        void set_cur_time (double cur_time);
        void set_cell_cycle_start_time ();
        void check_if_replication ();
        void check_if_division (RNG &generator);
        void replicate_genes ();
        void perform_dosage_compensation ();
        void swap_rxn_order (rxn_order_struct &A, rxn_order_struct &B);
        void set_cell_cycle_frozen_state (bool is_frozen_cell_cycle = false,
            std::string frozen_state = "G2");
        void move_to_G1 (RNG &generator, bool init_dosage_comp_adj=true);
        void move_to_G2 (bool init_dosage_comp_adj=true);
        std::string get_cur_cell_cycle_state () override;
    };
}

#endif
