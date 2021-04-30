// Optimized gillespie header file
#ifndef SDM_H
#define SDM_H
/********************************************//**
 \brief A class for Gillespie's stochastic simulation algorithm.

 The Gillespie class creates an object to perform
 exact stochastic simulation for any given chemical
 reaction network.
 ***********************************************/

class Gillespie: public Graph {
private:
  /* data */
  /********************************************//**
  \brief Number of chemical reaction channels.

  This is the total number of unique chemical
  raction channels that are present in the  user
  provided chemical reaction network.
   ***********************************************/
  int num_rxns;

  /********************************************//**
  \brief Number of molecular species.

  This is the total number of unique molecular
  species that are present in the user provided
  chemical reaction network.
   ***********************************************/
  int num_species;

  /********************************************//**
  \brief Time period for simulation run.

  The total time period for which the chemical
  reaction network should be simulated.
   ***********************************************/
  double time;

  // Vector of current molecule counts for all species.
  /********************************************//**
   \brief Vector of current molecule counts for all species.

   This stores the current molecule counts for all
   species in the chemical reaction network along
   the simulation trajectory.
   ***********************************************/
  std::vector<int> molecule_count_cur;

  /********************************************//**
   \brief Save time series count data.

   This boolean variable decides whether time series
   count data should be stored or not. Default value
   is false, when only steady state count values are
   reported.
   ***********************************************/

  bool save_timeseries;

  // Number of previous time points for which molecular counts are stored.
  int num_timepoints_save;

  // Vector of vectors to store running molecule counts for all species for last 'num_timepoints_save' time points.
  std::vector<std::vector<int>> molecule_count_history;

  // Vector for time points along the simulation path.
  std::vector<double> time_history;

  // Vector of reaction rates for all rxn channels.
  std::vector<double> rxn_rates;

  // Data structure to store stoichiometry info. for all rxn channels.
  std::vector<map<int, map<int, std::vector<int>>>> rxn_stoi;

  // Vector of rxn propensities for all reaction channels.
  std::vector<double> rxn_propensity;

  // Vector for reaction firing frequency for all reaction channels.
  std::vector<int> rxn_freq;

public:
  /* Member functions */
  // Constructor
  Gillespie (int num_rxns, int num_species);

  // Virtual destructor
  virtual ~Gillespie ();

  // Add reaction rates
  void add_rxn_rate (int rxn_id, double rate);

  // Add rxn input stoichiometric information.
  void add_rxn_stoi_input (int rxn_id, int mol_id, int stoi);

  // Add rxn output stoichiometric information.
  void add_rxn_stoi_output (int rxn_id, int mol_id, int stoi);

  // Update molecule count.
  void update_molecule_count ();

  // Update rxn propensity.
  void update_propensity (int rxn_id);

  // Sample next time step.
  double sample_time_step ();

  // Sample next reaction id.
  int sample_next_rxn ();

  // Simulate.
  void simulate ();

  // Update reaction frequencies.
  void update_rxn_freq (int next_rxn);

  // Set time period for simulation run.
  void set_time (double time);

  // Initialize molecule counts.
  void initialize_molecule_count (int species_id, int count);

  // Set number of previous history counts to store.
  void set_num_history (bool save_timeseries, int num_timepoints_save);
};

#endif
