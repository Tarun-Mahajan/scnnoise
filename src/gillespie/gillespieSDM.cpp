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
      GillespieSDM::GillespieSDM (int num_rxns, int num_species, int num_nodes_GRN):
        GillespieSSA (num_rxns, num_species, num_nodes_GRN){
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
