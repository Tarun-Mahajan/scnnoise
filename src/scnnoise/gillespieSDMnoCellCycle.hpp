// Optimized gillespie header file
#ifndef SDM_H
#define SDM_H

#include <vector>
#include <string>
#include "gillespieSSA.hpp"

namespace ScnnoiseInterface {
  /********************************************//**
   \brief A class for Gillespie's stochastic simulation algorithm.

   The Gillespie class creates an object to perform
   exact stochastic simulation for any given chemical
   reaction network.
   ***********************************************/

  class GillespieSDM : public GillespieSSA {
  private:

  public:
    /* Member functions */
    // Constructor
    GillespieSDM (int num_rxns, int num_genes,
      std::vector<int> num_species_gene_type,
      std::vector<int> num_rxns_gene_type, double max_time,
      bool save_timeseries, int num_timepoints_save,
      std::string count_save_file);

    void sort_reaction (int &rxn_selected) override;
    void update_cell_cycle_state (double cur_time) override;
  };
}

#endif
