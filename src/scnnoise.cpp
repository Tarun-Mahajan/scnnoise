// scNNoiSE simulator class
#include <vector>
#include <algorithm>
#include <iostream>
#include "GRN_simulation/graph.hpp"
#include "scnnoise.hpp"


// Constructor
scNNoiSE::scNNoiSE (int num_rxns, int num_species, int num_nodes_GRN) {
  const int sz = 1000;
  this->num_rxns = num_rxns;
  this->num_species = num_species;
  this->num_nodes_GRN = num_nodes_GRN;
  molecule_count_cur.reserve(sz);
  rxn_order.reserve(sz);
  network.reserve(1);
  network[0] = graph(num_nodes_GRN);
  node_type.reserve(sz);
  num_splice_variants.reserve(sz);
  rxn_rates.reserve(sz);
  species_order_GRN.reserve(sz);

  for (int i = 0; i < num_rxns; ++i) {
    rxn_order.push_back(i);
  }
}


void scNNoiSE::init_molecular_count (int molecule_id, int molecule_count) {
  molecule_count_cur[molecule_id] = molecule_count;
}


void scNNoiSE::add_GRNnode_type (int node_id, std::string str_type) {
  node_type[node_id] = str_type;
}


void scNNoiSE::add_num_splice_variants (int gene_id, int num_AS) {
  num_splice_variants[gene_id] = num_AS;
}


void scNNoiSE::add_rxn_reactants (int rxn_id, std::vector<int> rxn_reactants,
std::vector<int> reactants_stoichio) {
  for (auto &rxn_react: rxn_reactants) {
    reactions[rxn_id]['reactants'].push_back(rxn_react);
  }

  for (auto &rxn_stoi: reactants_stoichio) {
    reactions[rxn_id]['reactant stoichiometry'].push_back(rxn_stoi);
  }
}

void scNNoiSE::add_rxn_products (int rxn_id, std::vector<int> rxn_products,
std::vector<int> products_stoichio) {
  for (auto &rxn_prod: rxn_products) {
    reactions[rxn_id]['products'].push_back(rxn_prod);
  }

  for (auto &rxn_stoi: products_stoichio) {
    reactions[rxn_id]['product stoichiometry'].push_back(rxn_stoi);
  }
}

void scNNoiSE::add_rxn_rates (int rxn_id, double rxn_rate) {
  rxn_rates[rxn_id] = rxn_rate;
}


void scNNoiSE::add_species_order_GRN (int molecule_id, int GRN_order) {
  species_order_GRN[molecule_id] = GRN_order;
}
