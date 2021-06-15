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

    class GillespieSSA: public scnnoise {
    private:
      /* data */

      /********************************************//**
      \brief Time period for simulation run.

      The total time period for which the chemical
      reaction network should be simulated.
       ***********************************************/
      double time;

      /********************************************//**
       \brief Save time series count data.

       This boolean variable decides whether time series
       count data should be stored or not. Default value
       is false, when only steady state count values are
       reported.
       ***********************************************/
      bool save_timeseries;

      /********************************************//**
       \brief Number of previous time points for which molecular counts are stored in
              a vector of vectors..

       If save_timeseries is set to True, num_timepoints_save gives the number of
       past values of molecular counts to store in an 2D array
       (vector of vectors) of size . Once num_timepoints_save past values have been
       stored for each moleulcar species, the contents of the array are
       appened to an output file, and the array is emptied. After this,
       array starts storing the futures values of molecular counts until
       num_timepoints_save are stored, when the values are again appened to the
       output file. This process continues until the simulation terminates.
       At termination, the array again appends its values to the output file.
       ***********************************************/
      int num_timepoints_save;

      /********************************************//**
       \brief Vector of vectors to store running molecule counts for all species
              for last 'num_timepoints_save' time points.

       A 2D array stored as a vector of vectors to store running molecular counts
       for all the molecular specis in the system. This array ensures that the
       molecular counts are written to an output file at a fixed interval. If
       num_timepoints_save > 1, This reduces the number of times data has to be
       written to the output file.
       ***********************************************/
      std::vector<std::vector<int>> molecule_count_history;

      /********************************************//**
       \brief Vector to store time points along the simulation path.
       ***********************************************/
      std::vector<double> time_history;

      int num_threads;

      /********************************************//**
       \brief Vector to store reaction propensities.

       If \f$k\f$ is the reaction rate for a chemical reaction, then propensity,
       \f$P = kN\f$, where \f$N\f$ is the number of different ways reactant species
       can interact. \f$N\f$ is a function of the stoichiometric coefficients of the
       reactants.
       ***********************************************/
      std::vector<double> rxn_propensity;

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
