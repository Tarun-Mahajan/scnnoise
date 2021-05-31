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

/* function definitions */
// Constructor
GillespieSDM::GillespieSDM (int num_rxns, int num_species, int num_nodes_GRN):
  scnnoise (num_rxns, num_species, num_nodes_GRN){
  const int sz = 1000;
  this->num_rxns = num_rxns;
  this->num_species = num_species;
  time = 0;
  molecule_count_cur.reserve(sz);
  save_timeseries = false;
  num_timepoints_save = 1000;
  time_history.reserve(int(5*sz));
  time_history.push_back(0);
  rxn_rates.reserve(sz);
  rxn_stoi.reserve(sz);
  rxn_propensity.reserve(sz);
  rxn_freq.reserve(sz);
}

inline double GillespieSDM::sample_time_step (RNG &generator) {
  thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
  double rand_num = distribution(generator);
  return -log(1 - rand_num)/total_propensity;
}

int GillespieSDM::sample_next_rxn (RNG &generator) {
  thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
  double rand_num = distribution(generator);
  double selector = total_propensity * (double(1.0) - rand_num);
  int rxn_selected = -1;

  for ( int i = 0; i < num_rxns; ++i) {
    selector -= rxn_order[i].propensity_val;
    if (selector <= 0) {
      rxn_selected = i;
      break;
    }
  }
  if (rxn_selected > 0) {
    std::swap(rxn_order[rxn_selected - 1], rxn_order[rxn_selected]);
    rxn_selected -= 1;
  }
  return rxn_selected;
}

void GillespieSDM::update_fired_Reaction (int rxn_selected, bool &GRN_out_changed) {
  /*
  Find the gene and reaction type for the selected reaction channel. Extract the reactants,
  products and their stoichiometric coefficients for the selected reaction channel.
  */
  int gene_selected = rxn_order[rxn_selected].gene_id;
  // int rxn_type_selected = rxn_order[rxn_selected].rxn_type;
  int rxn_index = rxn_order[rxn_selected].rxn_type;
  // std::vector<int> rxn_temp = reactions[gene_selected].rxn_type;
  // std::vector<int>::iterator it = std::find(rxn_temp.begin(), rxn_temp.end(),
  //   rxn_type_selected);
  // int rxn_index = std::distance(rxn_temp.begin(), it);
  std::vector<int> rxn_selected_reactants =
    reactions[gene_selected].rxns[rxn_index].reactants;
  std::vector<int> rxn_selected_products =
    reactions[gene_selected].rxns[rxn_index].products;
  std::vector<int> rxn_selected_reactants_stoichio =
    reactions[gene_selected].rxns[rxn_index].reactants_stoichio;
  std::vector<int> rxn_selected_products_stoichio =
    reactions[gene_selected].rxns[rxn_index].products_stoichio;

  /*
  Subtract current reaction propensity from total reaction propensity for
  computing updated total propensity later.
  */
  total_propensity -= reactions[gene_selected].rxns[rxn_index].propensity_val;

  /*
  Update molecular count for reactants which change by firing the selected
  reaction channel.
  */
  int count_not_changed_reactants = 0;
  vector<bool> flag_changed_product_count(rxn_selected_products.size(), false);
  for (auto r = rxn_selected_reactants.begin(); r != rxn_selected_reactants.end(); ++r) {
    int reactant_index = std::distance(rxn_selected_reactants.begin(), r);
    std::vector<int>::iterator it =
    std::find(rxn_selected_products.begin(), rxn_selected_products.end(),
    *r);
    if (it != rxn_selected_products.end()) {
      int product_index = std::distance(rxn_selected_products.begin(), it);
      flag_changed_product_count[product_index] = true;
      int count_change = rxn_selected_products_stoichio[product_index] -
        rxn_selected_reactants_stoichio[reactant_index]
      if (count_change) {
        reactions[gene_selected].molecule_count_cur[*r] += count_change;
        if (*r == reactions[gene_selected].GRN_rxn_OUT) {
          GRN_out_changed = true;
        }
      }else{
        count_not_changed_reactants += 1;
      }
    }else{
      reactions[gene_selected].molecule_count_cur[*r] -=
        rxn_selected_reactants_stoichio[reactant_index];
    }
  }

  /*
  Update molecular count for products which change by firing the selected
  reaction channel.
  */
  for (auto r = rxn_selected_products.begin(); r != rxn_selected_products.end(); ++r) {
    int product_index = std::distance(rxn_selected_products.begin(), r);
    if (!flag_changed_product_count[product_index]) {
      reactions[gene_selected].molecule_count_cur[*r] +=
        rxn_selected_products_stoichio[product_index];
        if (*r == reactions[gene_selected].GRN_rxn_OUT) {
          GRN_out_changed = true;
        }
    }
  }

  /*
  Update propensity for the selected reaction channel.
  */
  if (count_not_changed_reactants < rxn_selected_reactants.size()) {
    double new_propensity = reactions[gene_selected].rxns[rxn_index].rxn_rate;
    for (auto r = rxn_selected_reactants.begin(); r != rxn_selected_reactants.end(); ++r) {
      int reactant_index = std::distance(rxn_selected_reactants.begin(), r);
      int reactant_stoichio = rxn_selected_reactants_stoichio[reactant_index];
      new_propensity *= (factorial(reactions[gene_selected].molecule_count_cur[*r])/
        factorial(reactions[gene_selected].molecule_count_cur[*r] - reactant_stoichio));
    }
    if (rxn_index == reactions[gene_selected].GRN_rxn_IN) {
      new_propensity *= regulation_function(gene_selected);
    }
    reactions[gene_selected].rxns[rxn_index].propensity_val = new_propensity;
    rxn_order[rxn_selected].propensity_val = new_propensity;

  }
  total_propensity += reactions[gene_selected].rxns[rxn_index].propensity_val;
}

void GillespieSDM::update_dependent_count_propensity (int rxn_selected,
  bool &GRN_out_changed) {

    /*
    Update propensity for dependent reactions belonging to the same gene as the
    reaction selected for firing
    */
    int gene_selected = rxn_order[rxn_selected].gene_id;
    int rxn_index = rxn_order[rxn_selected].rxn_type;
    int gene_type = reactions[gene_selected].gene_type;
    vector<int> rxn_selected_children;
    rxn_selected_children.reserve(gene_rxn_dependency[gene_type].get_size());
    gene_rxn_dependency[gene_type].find_children(rxn_index, rxn_selected_children);
    std::vector<int> rxn_selected_reactants;
    rxn_selected_reactants.reserve(reactions[gene_selected].molecule_count_cur.size());
    std::vector<int> rxn_selected_reactants_stoichio;
    rxn_selected_reactants_stoichio.reserve(reactions[gene_selected].molecule_count_cur.size());

    for (auto &rxn : rxn_selected_children) {
      total_propensity -= reactions[gene_selected].rxns[rxn].propensity_val;
      rxn_selected_reactants =
        reactions[gene_selected].rxns[rxn].reactants;
      rxn_selected_reactants_stoichio =
        reactions[gene_selected].rxns[rxn].reactants_stoichio;

      double new_propensity = reactions[gene_selected].rxns[rxn].rxn_rate;
      for (auto r = rxn_selected_reactants.begin(); r != rxn_selected_reactants.end(); ++r) {
        int reactant_index = std::distance(rxn_selected_reactants.begin(), r);
        int reactant_stoichio = rxn_selected_reactants_stoichio[reactant_index];
        new_propensity *= (factorial(reactions[gene_selected].molecule_count_cur[*r])/
          factorial(reactions[gene_selected].molecule_count_cur[*r] - reactant_stoichio));
      }
      if (rxn == reactions[gene_selected].GRN_rxn_IN) {
        new_propensity *= regulation_function(gene_selected);
      }
      reactions[gene_selected].rxns[rxn].propensity_val = new_propensity;
      for (auto &r : rxn_order) {
        if ((r.gene_id == gene_selected) && (r.rxn_type == rxn)) {
          r.propensity_val = new_propensity;
          break;
        }
      }
      total_propensity += reactions[gene_selected].rxns[rxn].propensity_val;

    }


    /*
    Update propensity for dependent reactions belonging to genes
    downstream of the gene related to the reaction selected for firing
    */
    if (GRN_out_changed) {
      vector<int> gene_children;
      network[0].find_children(gene_selected, gene_children);
      for (auto &g : gene_children) {
        int rxn = reactions[g].GRN_rxn_IN;
        total_propensity -= reactions[g].rxns[rxn].propensity_val;
        rxn_selected_reactants =
          reactions[g].rxns[rxn].reactants;
        rxn_selected_reactants_stoichio =
          reactions[g].rxns[rxn].reactants_stoichio;

        double new_propensity = reactions[g].rxns[rxn].rxn_rate;
        for (auto r = rxn_selected_reactants.begin(); r != rxn_selected_reactants.end(); ++r) {
          int reactant_index = std::distance(rxn_selected_reactants.begin(), r);
          int reactant_stoichio = rxn_selected_reactants_stoichio[reactant_index];
          new_propensity *= (factorial(reactions[g].molecule_count_cur[*r])/
            factorial(reactions[g].molecule_count_cur[*r] - reactant_stoichio));
        }
        new_propensity *= regulation_function(g);
        reactions[g].rxns[rxn].propensity_val = new_propensity;
        for (auto &r : rxn_order) {
          if ((r.gene_id == g) && (r.rxn_type == rxn)) {
            r.propensity_val = new_propensity;
            break;
          }
        }
        total_propensity += reactions[g].rxns[rxn].propensity_val;
      }
    }
}

double GillespieSDM::simulate (int num_threads) {
  compute_total_propensity();

  std::vector<std::uint_least32_t> rd_seeds = {rd(), rd(), rd(), rd()};
  std::seed_seq sd(rd_seeds.begin(), rd_seeds.end());
  thread_local static RNG generator{sd};
  bool stop_sim = false;

  while (!stop_sim) {
    double next_time_step = sample_time_step(generator);
    int next_rxn = sample_next_rxn(generator);
    bool GRN_out_changed = false;
    double total_time = std::accumulate(time_history.begin(), time_history.end(),
      decltype(time_history)::value_type(0)) + next_time_step;
    if (total_time < max_time) {
      update_fired_Reaction(next_rxn, GRN_out_changed);
      update_dependent_count_propensity(next_rxn, GRN_out_changed);
      time_history.push_back(next_time_step);

    }else{
      stop_sim = true;
    }
  }
}
