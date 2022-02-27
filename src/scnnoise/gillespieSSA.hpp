// Gillespie SSA interface header file
#ifndef SSA_H
#define SSA_H

#include <vector>
#include <random>
#include "scnnoise.hpp"

namespace ScnnoiseInterface {
    /********************************************//**
    \brief A class for Gillespie's stochastic simulation algorithm.

    The Gillespie class creates an object to perform
    exact stochastic simulation for any given chemical
    reaction network.
    ***********************************************/

    class GillespieSSA: public scNNoiSE {
    private:
    /* data */

    public:
        // typedef std::mt19937 RNG;
        /* Member functions */
        // Constructor
        GillespieSSA (int num_genes, std::string gene_filepath,
            std::string molecule_count_filepath,
            std::string count_save_file, bool keep_GRN,
            std::string GRN_filepath, int num_timepoints_save);

        // Sample next time step.
        double sample_time_step (RNG &generator);

        // Sample next reaction id.
        int sample_next_rxn (RNG &generator);

        // Update molecule counts by firing the selected reaction
        std::vector<std::string> update_fired_reaction (int rxn_selected);

        // Update propensity for reactions dependent via the GRN
        void update_dependent_count_propensity (int rxn_selected,
            std::vector<std::string> &GRN_out_changed);

        // Simulate.
        void simulate (bool compute_statistics = false,
            std::string statistics_file = "dummy", bool verbose = true) override;

        virtual void sort_reaction (int &rxn_selected) = 0;

        virtual void update_cell_cycle_state (double next_time,
            double cur_time, RNG &generator) = 0;

        virtual std::string get_cur_cell_cycle_state () = 0;

        virtual void init_cell_cycle_state (RNG &generator, double cur_time) = 0;

        // Update molecule count for reactants belonging to the fired reaction
        // channel
        // void update_fired_reaction_reactants(int gene_selected,
        //   int rxn_index, int &count_not_changed_reactants,
        //   std::vector<bool> &flag_changed_product_count,
        //   std::vector<bool> &GRN_out_changed,
        //   std::vector<int> &rxn_selected_reactants,
        //   std::vector<int> &rxn_selected_products,
        //   std::vector<int> &rxn_selected_reactants_stoichio,
        //   std::vector<int> &rxn_selected_products_stoichio);

        // Write column names for count_save_file
        void start_molecule_count_history_file ();

        // Update molecule count history
        void update_molecule_count_history (int &num_history, int &num_save_loop,
            bool simulation_ended, double cur_time);

        void save_molecule_count_at_interval (double time_prev, double time_next,
            double &points_collected_prev);

        void start_statistics_file (std::string statistics_file);

        void set_size_statistics_containers ();

        void set_size_statistics_history_containers ();

        void upate_running_statistics (double total_time_prev, double next_time_step);

        void write_statistics_to_file (std::string statistics_file);

        void save_molecule_count_at_random_times (double time_prev, double time_next,
            unsigned int &which_random_time_saved);
    };
}

#endif
