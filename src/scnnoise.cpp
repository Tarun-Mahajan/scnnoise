// scNNoiSE simulator class
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include "GRN_simulation/graph.hpp"
#include "scnnoise.hpp"


namespace ScnnoiseInterface {
  // Constructor
  scNNoiSE::scNNoiSE (int num_rxns, int num_species, int num_genes,
    const std::vector<int> num_species_gene_type,
    const std::vector<int> num_rxns_gene_type) {
    this->num_rxns = num_rxns;
    this->num_species = num_species;
    this->num_genes = num_genes;
    rxn_order.reserve(num_rxns);
    network.reserve(1);
    network[0] = graph(num_genes);
    reactions.reserve(num_genes);

    // num_species_gene_type.resize(num_gene_types);
    // num_species_gene_type.assign({4, 5, 4, 7});
    //
    // num_rxn_gene_type.resize(num_gene_types);
    // num_rxn_gene_type.assign({6, 7, 5, 9});

    this->num_species_gene_type.reserve(sz);
    this->num_species_gene_type = num_species_gene_type;
    this->num_rxns_gene_type.reserve(num_rxns_gene_type.size());
    this->num_rxns_gene_type = num_rxns_gene_type;

    /********************************************//**
     \brief Initialize dependency graph for different gene types.
     ***********************************************/
    gene_rxn_dependency.reserve(num_rxns_gene_type.size());
    int count = 0;
    for (auto &n : num_rxns_gene_type) {
      gene_rxn_dependency.push_back(graph(n));
    }
  }


  void scNNoiSE::init_molecular_count (int molecule_id, int molecule_count) {
    molecule_count_cur[molecule_id] = molecule_count;
  }

  void scNNoiSE::add_gene_state (int gene_id, int gene_type, int num_splice_variants) {
      reactions[gene_id].gene_type = gene_type;
      reactions[gene_id].num_splice_variants = num_splice_variants;
    }

  void scNNoiSE::add_reaction ()

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

  void scNNoiSE::add_GRN_edge (int src, int dest) {
    network[0].add_edge(src, dest);
    network[0].add_parent(src, dest);
  }

  void scNNoiSE::add_dependency_edge (int gene_type, int src, int dest) {
    gene_rxn_dependency[gene_type][0].add_edge(src, dest);
  }

  int scNNoiSE::factorial (int num) {
    if ((num == 0) || (num == 1)) {
      return 1;
    }else {
      return num*factorial(num - 1);
    }
  }

  void scNNoiSE::compute_total_propensity () {
    for (auto &rxn : rxn_order) {
      total_propensity += rxn.propensity_val;
    }
  }
}
