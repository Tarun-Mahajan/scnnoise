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
        {{"alpha", 12.0}, {"beta", 1.0/0.0083}};
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
            std::string GRN_filepath);

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

        void init_cell_cycle_state (RNG &generator, double cur_time) override;
        void cell_division (RNG& generator);
        void remove_dosage_compensation ();
        void update_propensity_cell_cycle ();
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
        void move_to_G1 (RNG &generator);
        void move_to_G2 ();
    };
}

#endif
