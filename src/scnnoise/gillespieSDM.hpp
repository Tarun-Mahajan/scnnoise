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
        GillespieSDM (int num_rxns, int num_genes,
          const std::vector<int> num_species_gene_type,
          const std::vector<int> num_rxns_gene_type, double max_time,
          bool save_timeseries, int num_timepoints_save);

        // Virtual destructor
        virtual ~GillespieSDM ();

        void sort_reaction (int &rxn_selected);
      };
    }
  }
}

#endif
