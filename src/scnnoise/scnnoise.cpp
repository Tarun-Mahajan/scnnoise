// scNNoiSE simulator class
#include "graph.hpp"
#include "scnnoise.hpp"
#include "graph_derived.hpp"
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <math.h>
// #include "graph.hpp"



namespace ScnnoiseInterface {
  // Constructor
  scNNoiSE::scNNoiSE (int num_rxns, int num_genes,
                     std::vector<int> num_species_gene_type,
                     std::vector<int> num_rxns_gene_type, double max_time,
                     bool save_timeseries, int num_timepoints_save, std::string count_save_file) {
                    this->num_rxns = num_rxns;
                    this->num_genes = num_genes;
                    rxn_order.reserve(num_rxns);
                    network.reserve(1);
                    network.push_back(GraphSpace::GRN(num_genes));
                    reactions.reserve(num_genes);
                    total_propensity = 0;

                    // num_species_gene_type.resize(num_gene_types);
                    // num_species_gene_type.assign({4, 5, 4, 7});
                    //
                    // num_rxn_gene_type.resize(num_gene_types);
                    // num_rxn_gene_type.assign({6, 7, 5, 9});
                    this->max_time = max_time;
                    this->save_timeseries = save_timeseries;
                    this->num_timepoints_save = num_timepoints_save;
                    this->num_species_gene_type.reserve(num_species_gene_type.size());
                    this->num_species_gene_type = num_species_gene_type;
                    this->num_rxns_gene_type.reserve(num_rxns_gene_type.size());
                    this->num_rxns_gene_type = num_rxns_gene_type;
                    this->count_save_file = count_save_file;

                    /********************************************//**
                     \brief Initialize dependency graph for different gene types.
                     ***********************************************/
                    gene_rxn_dependency.reserve(num_rxns_gene_type.size());
                    int count = 0;
                    for (auto &n : num_rxns_gene_type) {
                      gene_rxn_dependency.push_back(GraphSpace::GraphDependency(n));
                    }


                    molecule_count_history.reserve(num_genes);
                    for (int gene; gene < num_genes; ++gene) {
                      std::vector<std::vector<int>> count_gene;
                      count_gene.reserve(num_species_gene_type[gene]);
                      for (int id; id < num_species_gene_type[gene]; ++id) {
                        std::vector<int> count_species(num_timepoints_save, 0);
                        count_gene.push_back(count_species);
                      }
                      molecule_count_history.push_back(count_gene);
                    }

                    time_history.push_back(0);
  }


  void scNNoiSE::add_gene_state (int gene_id, int gene_type, std::vector<int> GRN_rxn_IN,
    std::vector<int> GRN_species_OUT, std::vector<int> molecule_count_cur,
    std::vector<std::vector<int>> reactants, std::vector<std::vector<int>> products,
    std::vector<std::vector<int>> reactants_stoichio, std::vector<std::vector<int>> products_stoichio,
    std::vector<double> rxn_rate, std::vector<double> propensity_val) {
      gene_rxn_channel_struct gene_state; // struct for gene state info
      std::vector<rxn_struct> rxns(num_rxns_gene_type[gene_id]); // struct for all rxns for a gene
      gene_state.gene_type = gene_type;
      gene_state.GRN_rxn_IN = GRN_rxn_IN;
      gene_state.GRN_species_OUT = GRN_species_OUT;
      gene_state.molecule_count_cur = molecule_count_cur;
      int rxn_index = 0;
      for (auto &rxn: rxns) {
        rxn.reactants = reactants[rxn_index];
        rxn.products = products[rxn_index];
        rxn.reactants_stoichio = reactants_stoichio[rxn_index];
        rxn.products_stoichio = products_stoichio[rxn_index];
        rxn.rxn_rate = rxn_rate[rxn_index];
        rxn.propensity_val = propensity_val[rxn_index];
        rxn_index += 1;
      }
      gene_state.rxns = rxns;
      reactions.push_back(gene_state);

      for (int id = 0; id < num_species_gene_type[gene_id]; ++id) {
        molecule_count_history[gene_id][id][0] = molecule_count_cur[id];
      }
  }


  void scNNoiSE::add_GRN_edge (int src, int dest, double prob_contr,
    double hill_coeff, double half_maximal, int rxn_IN, int species_OUT,
                              bool activator) {
                                  network[0].add_edge_kinetics(src, dest,
                                                              prob_contr,
                                                              hill_coeff,
                                                              half_maximal,
                                                              rxn_IN, species_OUT,
                                                              activator);
  }

  void scNNoiSE::add_dependency_edge (int gene_type, int src, int dest) {
    gene_rxn_dependency[gene_type].add_edge(src, dest);
  }

  int scNNoiSE::factorial (int num) {
    if ((num == 0) || (num == 1)) {
      return 1;
    }else {
      return num*factorial(num - 1);
    }
  }

  int scNNoiSE::factorial_ratio_propensity_func (int N, int r) {
      int propensity_factor = 1;
      for (int i = 0; i < r; ++i) {
          propensity_factor *= propensity_factor * (N - i);
      }
      return propensity_factor;
  }

  void scNNoiSE::compute_total_propensity () {
    for (auto &rxn : rxn_order) {
      total_propensity += rxn.propensity_val;
    }
  }

  double scNNoiSE::hill_function (int tf_count, double hill_coeff, double half_maximal,
                       bool activator) {
                         double tf_count_pow = pow(double(tf_count), hill_coeff);
                         double half_maximal_pow = pow(double(half_maximal), hill_coeff);
                         if (activator) {
                           return tf_count_pow/(tf_count_pow + half_maximal_pow);
                         }else{
                           return half_maximal_pow/(tf_count_pow + half_maximal_pow);
                         }
  }

  double scNNoiSE::regulation_function (int gene_selected, int rxn) {
    double regulation_val = 0;
    for (int src = 0; src < network[0].parent_list[gene_selected].size(); ++src) {
      std::vector<int>::iterator it =
        std::find(network[0].adj_list[src].begin(),
        network[0].adj_list[src].end(),
        gene_selected);
      int gene_selected_id = std::distance(network[0].adj_list[src].begin(), gene_selected);
      if (network[0].edge_rxn_params[src][gene_selected_id].rxn_IN == rxn) {
        int out_species = network[0].edge_rxn_params[src][gene_selected_id].species_OUT;
        double hill_coeff = network[0].edge_rxn_params[src][gene_selected_id].hill_coeff;
        double half_maximal = network[0].edge_rxn_params[src][gene_selected_id].half_maximal;
        double prob_contr = network[0].edge_rxn_params[src][gene_selected_id].prob_contr;
        bool activator = network[0].edge_rxn_params[src][gene_selected_id].activator;
        int tf_count = reactions[network[0].parent_list[gene_selected][src]].molecule_count_cur[out_species];
        regulation_val += prob_contr * hill_function(tf_count, hill_coeff, half_maximal,
                                                    activator);
      }
    }
    return regulation_val;
  }

  void init_gene_type_info () {
      /********************************************//**
       \brief  Classic two-state model of gene expression.
               No nascent and mature mRNA

       rxn_ids:
       0 - promoter in off state
       1 - promoter in on state
       2 - mRNA
       ***********************************************/
      std::vector<gene_type_struct> two_state(1);
      two_state
  }
  void swap_rxn_rates (std::vector<std::vector<int>> rxn_rates){
      //mbe add exception handling for rxn_rates size
      for(int gene = 0; gene<rxn_rates.size(); ++gene){
          for(int rxn = 0; rxn<rxn_rates[gene].size(); ++rxn){
            reactions[gene].rxns[rxn].rxn_rate = rxn_rates[gene][rxn];
          }
      }
  }

}
