#include <vector>
#include <algorithm>
#include <iostream>
#include "../GRN_simulation/graph.hpp"
#include "gillespie.hpp"

/* function definitions */
// Constructor
Gillespie::Gillespie (int num_rxns, int num_species):: Graph (num_rxns){
  this->num_rxns = num_rxns;
  this->num_species = num_species;
  time = 0;
  molecule_count_cur.resize(num_species);
  save_timeseries = false;
  num_timepoints_save = 100;
  time_history.push_back(0);
  rxn_rates.resize(num_rxns);
  rxn_stoi.resize(num_rxns);
  rxn_propensity.resize(num_rxns);
  rxn_freq.resize(num_rxns);
}

Gillespie::~Gillespie () {

}
