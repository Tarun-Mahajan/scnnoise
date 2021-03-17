// Optimized gillespie header file
#ifndef ODM_H
#define ODM_H
class Gillespie: public Graph {
private:
  /* data */
  // Number of reaction channels.
  int num_rxns;

  // Number of molecular species.
  int num_species;

  // Time period for simulation run.
  double time;

  // Vector of current molecule counts for all species.
  std::vector<int> molecule_count_cur;

  // Number of previous time points for which molecular counts are stored.
  int num_timepoints_save;

  // Vector of vectors for running molecule counts for all species for last 'num_timepoints_save' time points.
  std::vector<std::vector<int>> molecule_count_history;

  // Vector for time points along the simulation path.
  std::vector<double> time_history;

  // Vector of reaction rates for all rxn channels.
  std::vector<int> rxn_rates;

  // Data structure to store stoichiometry info. for all rxn channels.
  std::vector<map<int, map<int, std::vector<double>>>> rxn_stoi;

  // Vector of rxn propensities.
  std::vector<double> rxn_propensity;

  // Vector for reaction firing frequency.
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
  void update_molecule_count();

  // Update rxn propensity.
  void update_propensity (int rxn_id);

  // Sample next time step.
  double sample_time_step ();

  // Sample next reaction id.
  int sample_next_rxn();

  // Simulate.
  void simulate();

  // Update reaction frequencies.
  void update_rxn_freq(int next_rxn);
};

#endif
