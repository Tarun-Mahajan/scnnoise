// scNNoiSE simulator class
#include "graph.hpp"
#include "scnnoise.hpp"
#include "graph_derived.hpp"
#include "utils.hpp"
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

                    // Create dependency graphs
                    create_init_gene_type_info();

                    // Initialize gene states from file
                    init_gene_states_from_file (gene_filepath);

                    // Create GRN from input file
                    create_GRN(GRN_filepath);
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

  std::string scNNoiSE::match_and_return_gene_type (std::string in_gene_type) {
      std::string gene_return_type;
      bool is_found = false;
      for (auto const & it : gene_type_info) {
          in_gene_type_lower = in_gene_type;
          std::transform(in_gene_type_lower.begin(), in_gene_type_lower.end(),
            in_gene_type_lower.begin(), ::tolower);
          in_gene_type_lower_nospace =
            in_gene_type_lower.erase(std::remove_if(in_gene_type_lower.begin(),
            in_gene_type_lower.end(), ::isspace), in_gene_type_lower.end());

          gene_type = it.first;
          gene_type_lower = gene_type;
          std::transform(gene_type_lower.begin(), gene_type_lower.end(),
            gene_type_lower.begin(), ::tolower);
          gene_type_lower_nospace =
            gene_type_lower.erase(std::remove_if(gene_type_lower.begin(),
            gene_type_lower.end(), ::isspace), gene_type_lower.end());
          if ((gene_type == in_gene_type) ||
              (gene_type_lower == in_gene_type_lower) ||
              (gene_type_lower_nospace == in_gene_type_lower_nospace)) {
              is_found = true;
              gene_return_type = it.first;
              break;
          }
      }
      return gene_return_type;
  }

  void scNNoiSE::init_gene_states_from_file (std::string filepath) {
      std::ifstream gene_state_file(filepath);
      std::string row_text;
      std::string gene_name;
      std::string gene_type;
      std::vector<std::string>  GRN_rxn_IN;
      std::vector<std::string>  GRN_species_OUT;
      std::vector<std::string> rxn_names;
      std::vector<int> rxn_rates;
      std::string word;
      unsigned int num_rxn_IN = 1;
      unsigned int num_species_OUT = 0;
      unsigned int gene_count = 0;
      while (std::getline(gene_state_file, row_text)) {
          gene_rxn_channel_struct gene_rxns;
          GRN_rxn_IN.clear();
          GRN_species_OUT.clear();
          rxn_names.clear();
          rxn_rates.clear();
          GRN_activation.clear();
          std::istringstream str_stream(row_text);

          unsigned int id_counter = 0;
          while (std::getline(str_stream, word, ',')) {
              switch(id_counter) {
                  case 0:
                        {
                            gene_name = word;
                        }
                  case 1:
                        {
                            gene_type = word;
                        }
                  case 2:
                        {
                            num_rxn_IN = (unsigned int) std::stoi(word);
                        }
                  default:
                        {
                            if (id_counter > 2 && id_counter < 3 + num_rxn_IN) {
                                GRN_rxn_IN.push_back(word);
                            }else{
                                if (id_counter == 3 + num_rxn_IN) {
                                    num_species_OUT = (unsigned int) std::stoi(word);
                                }else{
                                    if (id_counter > 3 + num_rxn_IN &&
                                        id_counter <
                                        4 + num_rxn_IN + num_species_OUT) {
                                            GRN_species_OUT.push_back(word);
                                    }else{
                                        unsigned int diff_ = id_counter -
                                            (4 + num_rxn_IN + num_species_OUT);
                                        if (diff_ % 2 == 0) {
                                            rxn_names.push_back(word);
                                        }else{
                                            rxn_rates.push_back(std::stod(word));
                                        }
                                    }
                                }
                            }
                        }
              }
              ++id_counter;
          }
          gene_rxns.gene_name = gene_name;
          gene_map[gene_count] = gene_name;
          gene_rev_map[gene_name] = gene_count;
          gene_rxns.gene_type = match_and_return_gene_type(gene_type);
          gene_rxns.GRN_rxn_IN = GRN_rxn_IN;
          gene_rxns.GRN_species_OUT = GRN_species_OUT;
          for (unsigned int i = 0; i < rxn_names.size(); ++i) {
              gene_rxns.rxn_rates[rxn_names[i]] = rxn_rates[i];
          }
          gene_rxns.molecule_count_cur.resize(gene_type_info[gene_rxns.gene_type].num_species, 0);
          ++gene_count;
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

  void scNNoiSE::create_GRN (std::string filepath) {
      create_GRN_from_file(network[0], filepath, gene_rev_map);
  }

  // void scNNoiSE::add_dependency_edge (int gene_type, int src, int dest) {
  //   gene_rxn_dependency[gene_type].add_edge(src, dest);
  // }
  typedef std::map<std::string, std::map<std::string, int>> reactant_product_type;
  void scNNoiSE::add_new_dependency_graph (std::string gene_type_name,
    std::map<std::string, int> species_rev_map, std::map<int, std::string> rxn_map,
    reactant_product_type rxns_reactants,
    reactant_product_type rxns_products, std::vector<std::vector<int>> edge_list) {
        gene_type_struct new_gene_type;
        new_gene_type.species_rev_map = species_rev_map;
        for (auto const &it : species_rev_map) {
            new_gene_type.species_map[it.second] = it.first;
        }
        new_gene_type.num_species = species_rev_map.size();
        new_gene_type.rxn_map = rxn_map;
        for (auto const &it : rxn_map) {
            new_gene_type.rxn_rev_map[it.second] = it.first;
        }
        new_gene_type.num_rxns = rxn_map.size();
        for (auto const &it : rxns_reactants) {
            new_gene_type.rxns[it.first].reactants_stoichio = it.second;
        }
        for (auto const &it : rxns_products) {
            new_gene_type.rxns[it.first].products_stoichio = it.second;
        }
        new_gene_type.gene_rxn_dependency.push_back(
            GraphSpace::GraphDependency(new_gene_type.num_rxns));
        unsigned int node_count = 0;
        for (auto src : edge_list) {
            for (auto dest : src) {
                new_gene_type.gene_rxn_dependency[0].add_edge(node_count, dest);
            }
            ++node_count;
        }
        gene_type_info[gene_type_name] = new_gene_type;
    }

  gene_type_struct scNNoiSE::create_constitutive_type () {
      gene_type_struct gene_info;
      // Species are 0:gene, 1:mRNA, 2:protein
      gene_info.species_rev_map['gene'] = 0;
      gene_info.species_rev_map['mRNA'] = 1;
      gene_info.species_rev_map['protein'] = 2;
      for (auto const &it : gene_info.species_rev_map) {
          gene_info.species_map[it.second] = it.first;
      }
      gene_info.num_species = gene_info.species_map.size();
      // Reactions are 0:transcription, 1:mRNA decay, 2:translation,
      // 3:protein decay
      gene_info.rxn_map[0] = 'transcription';
      gene_info.rxn_map[1] = 'mRNA decay';
      gene_info.rxn_map[2] = 'translation';
      gene_info.rxn_map[3] = 'protein decay';
      for (auto const &it : gene_info.rxn_map) {
          gene_info.rxn_rev_map[it.second] = it.first;
      }
      gene_info.num_rxns = gene_info.rxn_map.size();
      std::map<std::string, rxn_struct> rxns_;
      // Transcription
      std::string str_ = "transcription";
      rxns_[str_].reactants_stoichio["gene"] = 1;
      rxns_[str_].products_stoichio["gene"] = 1;
      rxns_[str_].products_stoichio["mRNA"] = 1;
      // mRNA decay
      std::string str_ = "mRNA decay";
      rxns_[str_].reactants_stoichio["mRNA"] = 1;
      // Translation rxn
      std::string str_ = "translation";
      rxns_[str_].reactants_stoichio["mRNA"] = 1;
      rxns_[str_].products_stoichio["mRNA"] = 1;
      rxns_[str_].products_stoichio["protein"] = 1;
      // protein decay
      std::string str_ = "protein decay";
      rxns_[str_].reactants_stoichio["protein"] = 1;
      gene_info.rxns = rxns_;
      gene_info.gene_rxn_dependency.push_back(
          GraphSpace::GraphDependency(gene_info.num_rxns));
      gene_info.gene_rxn_dependency[0].add_edge(0, 1);
      gene_info.gene_rxn_dependency[0].add_edge(0, 2);
      gene_info.gene_rxn_dependency[0].add_edge(1, 1);
      gene_info.gene_rxn_dependency[0].add_edge(1, 2);
      gene_info.gene_rxn_dependency[0].add_edge(2, 3);
      gene_info.gene_rxn_dependency[0].add_edge(3, 3);
      return gene_info;
  }

  gene_type_struct scNNoiSE::create_constitutive_nascent_type () {
      gene_type_struct gene_info;
      // Species are 0:gene, 1:nascent mRNA, 2:mature mRNA, 3:protein
      gene_info.species_rev_map['gene'] = 0;
      gene_info.species_rev_map['nascent mRNA'] = 1;
      gene_info.species_rev_map['mature mRNA'] = 2;
      gene_info.species_rev_map['protein'] = 3;
      for (auto const &it : gene_info.species_rev_map) {
          gene_info.species_map[it.second] = it.first;
      }
      gene_info.num_species = gene_info.species_map.size();
      // Reactions are 0:transcription, 1:mRNA maturation, 2:mRNA decay, 3:translation,
      // 4:protein decay
      gene_info.rxn_map[0] = 'transcription';
      gene_info.rxn_map[1] = 'mRNA maturation';
      gene_info.rxn_map[2] = 'mRNA decay';
      gene_info.rxn_map[3] = 'translation';
      gene_info.rxn_map[4] = 'protein decay';
      for (auto const &it : gene_info.rxn_map) {
          gene_info.rxn_rev_map[it.second] = it.first;
      }
      gene_info.num_rxns = gene_info.rxn_map.size();
      std::map<std::string, rxn_struct> rxns_;
      // Transcription
      std::string str_ = "transcription";
      rxns_[str_].reactants_stoichio["gene"] = 1;
      rxns_[str_].products_stoichio["gene"] = 1;
      rxns_[str_].products_stoichio["nascent mRNA"] = 1;
      // mRNA maturation
      std::string str_ = "mRNA maturation";
      rxns_[str_].reactants_stoichio["nascent mRNA"] = 1;
      rxns_[str_].products_stoichio["mature mRNA"] = 1;
      // mRNA decay
      std::string str_ = "mRNA decay";
      rxns_[str_].reactants_stoichio["mature mRNA"] = 1;
      // Translation rxn
      std::string str_ = "translation";
      rxns_[str_].reactants_stoichio["mature mRNA"] = 1;
      rxns_[str_].products_stoichio["mature mRNA"] = 1;
      rxns_[str_].products_stoichio["protein"] = 1;
      // protein decay
      std::string str_ = "protein decay";
      rxns_[str_].reactants_stoichio["protein"] = 1;
      gene_info.rxns = rxns_;
      gene_info.gene_rxn_dependency.push_back(
          GraphSpace::GraphDependency(gene_info.num_rxns));
      gene_info.gene_rxn_dependency[0].add_edge(0, 1);
      gene_info.gene_rxn_dependency[0].add_edge(1, 1);
      gene_info.gene_rxn_dependency[0].add_edge(1, 2);
      gene_info.gene_rxn_dependency[0].add_edge(1, 3);
      gene_info.gene_rxn_dependency[0].add_edge(2, 2);
      gene_info.gene_rxn_dependency[0].add_edge(2, 3);
      gene_info.gene_rxn_dependency[0].add_edge(3, 4);
      gene_info.gene_rxn_dependency[0].add_edge(4, 4);
      return gene_info;
  }

  gene_type_struct scNNoiSE::create_two_state_type () {
      gene_type_struct gene_info;
      // Species are 0:gene off, 1:gene on, 2:mRNA, 3:protein
      gene_info.species_rev_map['gene off'] = 0;
      gene_info.species_rev_map['gene on'] = 1;
      gene_info.species_rev_map['mRNA'] = 2;
      gene_info.species_rev_map['protein'] = 3;
      for (auto const &it : gene_info.species_rev_map) {
          gene_info.species_map[it.second] = it.first;
      }
      gene_info.num_species = gene_info.species_map.size();
      // Reactions are 0:gene on, 1:gene off 2:transcription, 3:mRNA decay,
      // 4:translation, 5:protein decay
      gene_info.rxn_map[0] = 'gene on';
      gene_info.rxn_map[1] = 'gene off';
      gene_info.rxn_map[2] = 'transcription';
      gene_info.rxn_map[3] = 'mRNA decay';
      gene_info.rxn_map[4] = 'translation';
      gene_info.rxn_map[5] = 'protein decay';
      for (auto const &it : gene_info.rxn_map) {
          gene_info.rxn_rev_map[it.second] = it.first;
      }
      gene_info.num_rxns = gene_info.rxn_map.size();
      std::map<std::string, rxn_struct> rxns_;
      // Gene on
      std::string str_ = "gene on";
      rxns_[str_].reactants_stoichio["gene off"] = 1;
      rxns_[str_].products_stoichio["gene on"] = 1;
      // Gene off
      std::string str_ = "gene off";
      rxns_[str_].reactants_stoichio["gene on"] = 1;
      rxns_[str_].products_stoichio["gene off"] = 1;
      // Transcription
      std::string str_ = "transcription";
      rxns_[str_].reactants_stoichio["gene on"] = 1;
      rxns_[str_].products_stoichio["gene on"] = 1;
      rxns_[str_].products_stoichio["mRNA"] = 1;
      // mRNA decay
      std::string str_ = "mRNA decay";
      rxns_[str_].reactants_stoichio["mRNA"] = 1;
      // Translation rxn
      std::string str_ = "translation";
      rxns_[str_].reactants_stoichio["mRNA"] = 1;
      rxns_[str_].products_stoichio["mRNA"] = 1;
      rxns_[str_].products_stoichio["protein"] = 1;
      // protein decay
      std::string str_ = "protein decay";
      rxns_[str_].reactants_stoichio["protein"] = 1;
      gene_info.rxns = rxns_;
      gene_info.gene_rxn_dependency.push_back(
          GraphSpace::GraphDependency(gene_info.num_rxns));
      gene_info.gene_rxn_dependency[0].add_edge(0, 0);
      gene_info.gene_rxn_dependency[0].add_edge(0, 1);
      gene_info.gene_rxn_dependency[0].add_edge(0, 2);
      gene_info.gene_rxn_dependency[0].add_edge(1, 0);
      gene_info.gene_rxn_dependency[0].add_edge(1, 1);
      gene_info.gene_rxn_dependency[0].add_edge(1, 2);
      gene_info.gene_rxn_dependency[0].add_edge(2, 3);
      gene_info.gene_rxn_dependency[0].add_edge(2, 4);
      gene_info.gene_rxn_dependency[0].add_edge(3, 3);
      gene_info.gene_rxn_dependency[0].add_edge(3, 4);
      gene_info.gene_rxn_dependency[0].add_edge(4, 5);
      gene_info.gene_rxn_dependency[0].add_edge(5, 5);
      return gene_info;
  }

  gene_type_struct scNNoiSE::create_two_state_nascent_type () {
      gene_type_struct gene_info;
      // Species are 0:gene off, 1:gene on, 2:mRNA, 3:protein
      gene_info.species_rev_map['gene off'] = 0;
      gene_info.species_rev_map['gene on'] = 1;
      gene_info.species_rev_map['nascent mRNA'] = 2;
      gene_info.species_rev_map['mature mRNA'] = 3;
      gene_info.species_rev_map['protein'] = 4;
      for (auto const &it : gene_info.species_rev_map) {
          gene_info.species_map[it.second] = it.first;
      }
      gene_info.num_species = gene_info.species_map.size();
      // Reactions are 0:gene on, 1:gene off 2:transcription, 3:mRNA decay,
      // 4:translation, 5:protein decay
      gene_info.rxn_map[0] = 'gene on';
      gene_info.rxn_map[1] = 'gene off';
      gene_info.rxn_map[2] = 'transcription';
      gene_info.rxn_map[3] = 'mRNA maturation';
      gene_info.rxn_map[4] = 'mRNA decay';
      gene_info.rxn_map[5] = 'translation';
      gene_info.rxn_map[6] = 'protein decay';
      for (auto const &it : gene_info.rxn_map) {
          gene_info.rxn_rev_map[it.second] = it.first;
      }
      gene_info.num_rxns = gene_info.rxn_map.size();
      std::map<std::string, rxn_struct> rxns_;
      // Gene on
      std::string str_ = "gene on";
      rxns_[str_].reactants_stoichio["gene off"] = 1;
      rxns_[str_].products_stoichio["gene on"] = 1;
      // Gene off
      std::string str_ = "gene off";
      rxns_[str_].reactants_stoichio["gene on"] = 1;
      rxns_[str_].products_stoichio["gene off"] = 1;
      // Transcription
      std::string str_ = "transcription";
      rxns_[str_].reactants_stoichio["gene on"] = 1;
      rxns_[str_].products_stoichio["gene on"] = 1;
      rxns_[str_].products_stoichio["nascent mRNA"] = 1;
      // mRNA maturation
      std::string str_ = "mRNA maturation";
      rxns_[str_].reactants_stoichio["nascent mRNA"] = 1;
      rxns_[str_].products_stoichio["mature mRNA"] = 1;
      // mRNA decay
      std::string str_ = "mRNA decay";
      rxns_[str_].reactants_stoichio["mature mRNA"] = 1;
      // Translation rxn
      std::string str_ = "translation";
      rxns_[str_].reactants_stoichio["mature mRNA"] = 1;
      rxns_[str_].products_stoichio["mature mRNA"] = 1;
      rxns_[str_].products_stoichio["protein"] = 1;
      // protein decay
      std::string str_ = "protein decay";
      rxns_[str_].reactants_stoichio["protein"] = 1;
      gene_info.rxns = rxns_;
      gene_info.gene_rxn_dependency.push_back(
          GraphSpace::GraphDependency(gene_info.num_rxns));
      gene_info.gene_rxn_dependency[0].add_edge(0, 0);
      gene_info.gene_rxn_dependency[0].add_edge(0, 1);
      gene_info.gene_rxn_dependency[0].add_edge(0, 2);
      gene_info.gene_rxn_dependency[0].add_edge(1, 0);
      gene_info.gene_rxn_dependency[0].add_edge(1, 1);
      gene_info.gene_rxn_dependency[0].add_edge(1, 2);
      gene_info.gene_rxn_dependency[0].add_edge(2, 3);
      gene_info.gene_rxn_dependency[0].add_edge(3, 4);
      gene_info.gene_rxn_dependency[0].add_edge(3, 5);
      gene_info.gene_rxn_dependency[0].add_edge(4, 4);
      gene_info.gene_rxn_dependency[0].add_edge(4, 5);
      gene_info.gene_rxn_dependency[0].add_edge(5, 6);
      gene_info.gene_rxn_dependency[0].add_edge(6, 6);

      return gene_info;
  }

  void scNNoiSE::create_init_gene_type_info () {
      gene_type_info['constitutive'] = create_constitutive_type();
      gene_type_info['constitutive nascent'] = create_constitutive_type();
      gene_type_info['two-state'] = create_two_state_type();
      gene_type_info['two-state nascent'] = create_two_state_nascent_type();
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
}
