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
    typedef std::mt19937 RNG;
    /* Member functions */
    // Constructor
    GillespieSSA (int num_rxns, int num_genes,
      std::vector<int> num_species_gene_type,
      std::vector<int> num_rxns_gene_type, double max_time,
      bool save_timeseries, int num_timepoints_save,
      std::string count_save_file);

    // Sample next time step.
    double sample_time_step (RNG &generator);

    // Sample next reaction id.
    int sample_next_rxn (RNG &generator);

    // Update molecule counts by firing the selected reaction
    std::vector<bool> update_fired_reaction (int rxn_selected);

    // Update propensity for reactions dependent via the GRN
    void update_dependent_count_propensity (int rxn_selected,
                                           std::vector<std::string> &GRN_out_changed);

    // Simulate.
    void simulate () override;

    virtual void sort_reaction (int &rxn_selected) = 0;

    virtual void update_cell_cycle_state (double cur_time) = 0;

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
    void start_molecule_count_history_file ()

    // Update molecule count history
    void update_molecule_count_history (int &num_history, int &num_save_loop);
  };
}

#endif
