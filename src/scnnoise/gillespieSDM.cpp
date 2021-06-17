#include <vector>
#include <algorithm>
#include <iostream>
#include <omp.h>
#include <string>
#include <random>
#include <math.h>
#include "../GRN_simulation/graph.hpp"
#include "../scnnoise.hpp"
#include "gillespieSDM.hpp"
#include "gillespieSSA.hpp"

namespace ScnnoiseInterface {
  namespace SimulatorGillespieSSA {
    namespace SimulatorGillespieSDM {
      /* function definitions */
      // Constructor
      GillespieSDM::GillespieSDM (int num_rxns, int num_genes,
        const std::vector<int> num_species_gene_type,
        const std::vector<int> num_rxns_gene_type, double max_time,
        bool save_timeseries, int num_timepoints_save):
        GillespieSSA (num_rxns, num_genes, num_species_gene_type,
          num_rxns_gene_type, max_time, save_timeseries, num_timepoints_save) {
      }

      GillespieSDM::sort_reaction (int &rxn_selected) {
        if (rxn_selected > 0) {
          std::swap(rxn_order[rxn_selected - 1], rxn_order[rxn_selected]);
          rxn_selected -= 1;
        }
      }
    }
  }
}
