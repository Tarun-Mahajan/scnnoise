// Graph class and algorithms
#include "GRN.hpp"
#include "graph.hpp"
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>



namespace GraphSpace {
  // Constructor
  GRN::GRN (int N): Graph(N) {
    edge_rxn_params.reserve(N);
  }

  // Function to add edge
  void GRN::add_edge_kinetics (int src, int dest, double prob_contr,
                               double hill_coeff, double half_maximal,
                               int rxn_IN, int species_OUT) {
                               add_edge(src, dest);
                               edge_rxn_struct edge_rxn_temp;
                               edge_rxn_temp.prob_contr = prob_contr;
                               edge_rxn_temp.hill_coeff = hill_coeff;
                               edge_rxn_temp.half_maximal = half_maximal;
                               edge_rxn_temp.rxn_IN = rxn_IN;
                               edge_rxn_temp.species_OUT = species_OUT;
                               edge_rxn_params[src].push_back(edge_rxn_temp);
                             }
  inline void GRN::add_edge (int src, int dest) {
    adj_list[src].push_back(dest);
    parent_list[dest].push_back(src);
  }

  void GRN::find_children_edge_info (int vert,
    std::vector<edge_rxn_struct> &children_edge_info) {
      if (!edge_rxn_params[vert].empty()) {
        for (auto v: edge_rxn_params[vert]) {
          children_edge_info.push_back(v);
        }
      }
  }
}
