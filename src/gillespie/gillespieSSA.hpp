// Gillespie SSA interface header file
#ifndef SSA_H
#define SSA_H

namespace ScnnoiseInterface {
  namespace SimulatorGillespieSSA {
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
      GillespieSSA (int num_rxns, int num_species, int num_nodes_GRN, double time,
        bool save_timeseries, int num_timepoints_save
        std::vector<std::vector<int>> molecule_count_history,
        std::vector<double> time_history, std::vector<double> rxn_propensity);

      // Virtual destructor
      virtual ~GillespieSSA ();

      // Sample next time step.
      double sample_time_step ();

      // Sample next reaction id.
      int sample_next_rxn ();

      // Update molecule counts by firing the selected reaction
      std::vector<bool> update_fired_Reaction (int rxn_selected, bool &GRN_out_changed);

      // Update propensity for reactions dependent via the GRN
      void update_dependent_count_propensity (int rxn_selected, bool &GRN_out_changed);

      // Simulate.
      void simulate ();

      virtual void sort_reaction (int &rxn_selected) = 0;
    };
  }
}

#endif
