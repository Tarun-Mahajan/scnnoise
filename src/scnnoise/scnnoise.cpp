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
}

void scNNoiSE::init_molecular_count (int molecule_id, int molecule_count) {

}
