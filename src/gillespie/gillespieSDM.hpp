// Optimized gillespie header file
#ifndef SDM_H
#define SDM_H

namespace ScnnoiseInterface {
  namespace SimulatorGillespieSSA {
    namespace SimulatorGillespieSDM {
      /********************************************//**
       \brief A class for Gillespie's stochastic simulation algorithm.

       The Gillespie class creates an object to perform
       exact stochastic simulation for any given chemical
       reaction network.
       ***********************************************/

      class GillespieSDM: public GillespieSSA {
      private:

      public:
        typedef std::mt19937 RNG;
        /* Member functions */
        // Constructor
        GillespieSDM (int num_rxns, int num_species, int num_nodes_GRN, double time,
          bool save_timeseries, int num_timepoints_save
          std::vector<std::vector<int>> molecule_count_history,
          std::vector<double> time_history, std::vector<double> rxn_propensity);

        // Virtual destructor
        virtual ~GillespieSDM ();

        void sort_reaction (int &rxn_selected);
      };
    }
  }
}

#endif
