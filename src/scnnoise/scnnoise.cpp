// scNNoiSE simulator class
#include "graph.hpp"
#include "scnnoise.hpp"
#include "graph_derived.hpp"
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cctype>
#include <sstream>
#include <string>
#include <math.h>
// #include "graph.hpp"



namespace ScnnoiseInterface {
    // Constructor
    scNNoiSE::scNNoiSE (int num_genes, std::string gene_filepath,
        std::string GRN_filepath, std::string molecule_count_filepath,
        std::string count_save_file) {
        // this->num_rxns = num_rxns;
        this->num_genes = num_genes;
        // rxn_order.reserve(num_rxns);
        network.reserve(1);
        network.push_back(GraphSpace::GRN(num_genes));
        reactions.reserve(num_genes);
        total_propensity = 0;

        // num_species_gene_type.resize(num_gene_types);
        // num_species_gene_type.assign({4, 5, 4, 7});
        //
        // num_rxn_gene_type.resize(num_gene_types);
        // num_rxn_gene_type.assign({6, 7, 5, 9});
        this->max_time = 10000;
        this->save_timeseries = false;
        this->num_timepoints_save = 1000;
        // this->count_save_file = count_save_file;

        /********************************************//**
         \brief Initialize dependency graph for different gene types.
         ***********************************************/
        // gene_rxn_dependency.reserve(num_rxns_gene_type.size());
        // int count = 0;
        // for (auto &n : num_rxns_gene_type) {
        //   gene_rxn_dependency.push_back(GraphSpace::GraphDependency(n));
        // }


        time_history.push_back(0);

        // Create dependency graphs
        create_init_gene_type_info();

        // Initialize gene states from file
        init_gene_states_from_file(gene_filepath);

        // Create GRN from input file
        create_GRN(GRN_filepath);

        // Initialize max_rxn_rate_change
        init_max_rxn_rate_change();

        // Initialize molecule count
        init_molecule_count (molecule_count_filepath);

        // Initialize rxn order
        init_rxn_order();
    }

    void scNNoiSE::init_molecule_count (std::string filepath) {
        // for (auto const &it : molecule_count_init) {
        //     int gene_id = gene_rev_map[it.first];
        //     gene_rxn_channel_struct rxn_ = reactions[gene_id];
        //     gene_type_struct gene_info = gene_type_info[rxn_.gene_type];
        //     for (auto const &rxns : it.second) {
        //         reactions[gene_id].molecule_count_cur[gene_info.species_rev_map[rxns.first]] =
        //             rxns.second;
        //     }
        //
        // }


        std::ifstream count_state_file(filepath);
        std::string row_text;
        std::string gene_name;
        std::string gene_type;
        std::vector<std::string>  species_name;
        std::vector<int>  species_count;
        std::string word;
        while (std::getline(count_state_file, row_text)) {
            species_name.clear();
            species_count.clear();
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
                    default:
                        {
                            if ((id_counter - 2) % 2 == 0) {
                                species_name.push_back(word);
                            }else{
                                species_count.push_back(std::stoi(word));
                            }
                        }
                }
                ++id_counter;
            }
            int gene_id = gene_rev_map[gene_name];
            gene_rxn_channel_struct rxn_ = reactions[gene_id];
            gene_type_struct gene_info = gene_type_info[rxn_.gene_type];
            for (std::size_t i = 0; i < species_name.size(); ++i) {
                reactions[gene_id].molecule_count_cur[gene_info.species_rev_map[species_name[i]]] =
                    species_count[i];
            }

        }
    }

    void scNNoiSE::init_molecule_count_history () {
        molecule_count_history.resize(num_genes);
        for (int gene; gene < num_genes; ++gene) {
          std::vector<std::vector<int>> count_gene;
          count_gene.resize(reactions[gene].molecule_count_cur.size());
          for (std::size_t id = 0; id < reactions[gene].molecule_count_cur.size(); ++id) {
            std::vector<int> count_species(num_timepoints_save, 0);
            count_gene.push_back(count_species);
          }
          molecule_count_history.push_back(count_gene);
        }
    }

    std::string scNNoiSE::match_and_return_gene_type (std::string in_gene_type) {
        std::string gene_return_type;
        bool is_found = false;
        for (auto const & it : gene_type_info) {
            std::string in_gene_type_lower = in_gene_type;
            std::transform(in_gene_type_lower.begin(), in_gene_type_lower.end(),
                in_gene_type_lower.begin(), ::tolower);
            std::string in_gene_type_lower_nospace = in_gene_type_lower;
            in_gene_type_lower_nospace.erase(
                std::remove_if(in_gene_type_lower_nospace.begin(),
                in_gene_type_lower_nospace.end(), ::isspace),
                in_gene_type_lower_nospace.end());

            std::string gene_type = it.first;
            std::string gene_type_lower = gene_type;
            std::transform(gene_type_lower.begin(), gene_type_lower.end(),
                gene_type_lower.begin(), ::tolower);
                std::string gene_type_lower_nospace = gene_type_lower;
             gene_type_lower_nospace.erase(
                 std::remove_if(gene_type_lower_nospace.begin(),
                gene_type_lower_nospace.end(), ::isspace),
                gene_type_lower_nospace.end());
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

    void scNNoiSE::init_max_rxn_rate_change () {
        for (unsigned int i = 0; i < num_genes; ++i) {
            std::string gene_type = reactions[i].gene_type;
            for (auto const &it : gene_type_info[gene_type].rxn_rev_map) {
                max_rxn_rate_change[i][it.first] = 0.5;
            }

        }
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

    double scNNoiSE::compute_regulation_function (int gene_id, std::string rxn_name) {
        double regulation_function_factor = 0;
        gene_type_struct gene_info = gene_type_info[gene_map[gene_id]];
        int rxn_id = gene_info.rxn_rev_map[rxn_name];
        if (reactions[gene_id].GRN_rxn_IN.size() != 0){
            auto it = std::find(reactions[gene_id].GRN_rxn_IN.begin(),
                reactions[gene_id].GRN_rxn_IN.end(),
                rxn_name);

            if (it != reactions[gene_id].GRN_rxn_IN.end()) {
                std::vector<int> parents;
                network[0].find_parents(gene_id, parents);
                for (auto const &src : parents) {
                    for (std::size_t dest = 0; dest < network[0].adj_list[src].size(); ++dest) {
                        if ((network[0].adj_list[src][dest] == gene_id) &&
                            (network[0].edge_rxn_params[src][dest].rxn_IN == rxn_id)) {
                                int species_out_id = network[0].edge_rxn_params[src][dest].species_OUT;
                                int mol_count = reactions[src].molecule_count_cur[species_out_id];
                                regulation_function_factor +=
                                    hill_function(mol_count,
                                        network[0].edge_rxn_params[src][dest].hill_coeff,
                                        network[0].edge_rxn_params[src][dest].half_maximal,
                                        network[0].edge_rxn_params[src][dest].activator,
                                        network[0].edge_rxn_params[src][dest].prob_contr);
                        }
                    }
                }
            }
        }
        regulation_function_factor *= max_rxn_rate_change[gene_id][rxn_name];
        regulation_function_factor += 1;
        return regulation_function_factor;
    }

    double scNNoiSE::compute_propensity (std::string gene_name,
        std::string rxn_name) {
        int gene_id = gene_rev_map[gene_name];
        double new_propensity = reactions[gene_id].rxn_rates[rxn_name];
        std::map<std::string, int> rxn_stoi_map =
            gene_type_info[reactions[gene_id].gene_type].rxns[rxn_name].reactants_stoichio;
        if (rxn_stoi_map.size() > 0) {
            for (auto const &it : rxn_stoi_map) {
                int species_id = gene_type_info[reactions[gene_id].gene_type].species_rev_map[it.first];
                new_propensity *=
                    factorial_ratio_propensity_func(reactions[gene_id].molecule_count_cur[species_id],
                    it.second);
            }
        }

        new_propensity *= compute_regulation_function(gene_id, rxn_name);
        return new_propensity;
    }

    void scNNoiSE::init_rxn_order () {
        unsigned int order_count = 0;
        for (auto &it : reactions) {
            for (auto const &rxns_ : it.rxn_rates) {
                rxn_order_struct rxn_order_temp;
                rxn_order_temp.gene_id = gene_rev_map[it.gene_name];
                rxn_order_temp.rxn_name = rxns_.first;
                rxn_order_temp.propensity_val = compute_propensity(it.gene_name,
                    rxns_.first);
                it.propensity_vals[rxns_.first] = rxn_order_temp.propensity_val;
                rxn_order.push_back(rxn_order_temp);
                rxn_order_map[it.gene_name][rxns_.first] = order_count;
            }
        }
    }

    void scNNoiSE::create_GRN (std::string filepath) {
        this->create_GRN_from_file(filepath);
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
        for (auto &src : edge_list) {
            for (auto &dest : src) {
                new_gene_type.gene_rxn_dependency[0].add_edge(node_count, dest);
            }
            ++node_count;
        }
        gene_type_info[gene_type_name] = new_gene_type;
    }

    gene_type_struct scNNoiSE::create_constitutive_type () {
        gene_type_struct gene_info;
        // Species are 0:gene, 1:mRNA, 2:protein
        gene_info.species_rev_map["gene"] = 0;
        gene_info.species_rev_map["mRNA"] = 1;
        gene_info.species_rev_map["protein"] = 2;
        for (auto const &it : gene_info.species_rev_map) {
            gene_info.species_map[it.second] = it.first;
        }
        gene_info.num_species = gene_info.species_map.size();
        // Reactions are 0:transcription, 1:mRNA decay, 2:translation,
        // 3:protein decay
        gene_info.rxn_map[0] = "transcription";
        gene_info.rxn_map[1] = "mRNA decay";
        gene_info.rxn_map[2] = "translation";
        gene_info.rxn_map[3] = "protein decay";
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
        gene_info.species_rev_map["gene"] = 0;
        gene_info.species_rev_map["nascent mRNA"] = 1;
        gene_info.species_rev_map["mature mRNA"] = 2;
        gene_info.species_rev_map["protein"] = 3;
        for (auto const &it : gene_info.species_rev_map) {
            gene_info.species_map[it.second] = it.first;
        }
        gene_info.num_species = gene_info.species_map.size();
        // Reactions are 0:transcription, 1:mRNA maturation, 2:mRNA decay, 3:translation,
        // 4:protein decay
        gene_info.rxn_map[0] = "transcription";
        gene_info.rxn_map[1] = "mRNA maturation";
        gene_info.rxn_map[2] = "mRNA decay";
        gene_info.rxn_map[3] = "translation";
        gene_info.rxn_map[4] = "protein decay";
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
        gene_info.species_rev_map["gene off"] = 0;
        gene_info.species_rev_map["gene on"] = 1;
        gene_info.species_rev_map["mRNA"] = 2;
        gene_info.species_rev_map["protein"] = 3;
        for (auto const &it : gene_info.species_rev_map) {
            gene_info.species_map[it.second] = it.first;
        }
        gene_info.num_species = gene_info.species_map.size();
        // Reactions are 0:gene on, 1:gene off 2:transcription, 3:mRNA decay,
        // 4:translation, 5:protein decay
        gene_info.rxn_map[0] = "gene on";
        gene_info.rxn_map[1] = "gene off";
        gene_info.rxn_map[2] = "transcription";
        gene_info.rxn_map[3] = "mRNA decay";
        gene_info.rxn_map[4] = "translation";
        gene_info.rxn_map[5] = "protein decay";
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
        gene_info.species_rev_map["gene off"] = 0;
        gene_info.species_rev_map["gene on"] = 1;
        gene_info.species_rev_map["nascent mRNA"] = 2;
        gene_info.species_rev_map["mature mRNA"] = 3;
        gene_info.species_rev_map["protein"] = 4;
        for (auto const &it : gene_info.species_rev_map) {
            gene_info.species_map[it.second] = it.first;
        }
        gene_info.num_species = gene_info.species_map.size();
        // Reactions are 0:gene on, 1:gene off 2:transcription, 3:mRNA decay,
        // 4:translation, 5:protein decay
        gene_info.rxn_map[0] = "gene on";
        gene_info.rxn_map[1] = "gene off";
        gene_info.rxn_map[2] = "transcription";
        gene_info.rxn_map[3] = "mRNA maturation";
        gene_info.rxn_map[4] = "mRNA decay";
        gene_info.rxn_map[5] = "translation";
        gene_info.rxn_map[6] = "protein decay";
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

    gene_type_struct scNNoiSE::create_two_state_reduced_type () {
        gene_type_struct gene_info;
        // Species are 0:gene off, 1:gene on, 2:mRNA, 3:protein
        gene_info.species_rev_map["mRNA"] = 0;
        gene_info.species_rev_map["protein"] = 1;
        for (auto const &it : gene_info.species_rev_map) {
            gene_info.species_map[it.second] = it.first;
        }
        gene_info.num_species = gene_info.species_map.size();
        // Reactions are 0:gene on, 1:gene off 2:transcription, 3:mRNA decay,
        // 4:translation, 5:protein decay
        gene_info.rxn_map[0] = "transcription";
        gene_info.rxn_map[1] = "mRNA decay";
        gene_info.rxn_map[2] = "translation";
        gene_info.rxn_map[3] = "protein decay";
        for (auto const &it : gene_info.rxn_map) {
            gene_info.rxn_rev_map[it.second] = it.first;
        }
        gene_info.num_rxns = gene_info.rxn_map.size();
        std::map<std::string, rxn_struct> rxns_;
        // Transcription
        std::string str_ = "transcription";
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

    void scNNoiSE::create_init_gene_type_info () {
        gene_type_info["constitutive"] = create_constitutive_type();
        gene_type_info["constitutive nascent"] = create_constitutive_type();
        gene_type_info["two-state"] = create_two_state_type();
        gene_type_info["two-state nascent"] = create_two_state_nascent_type();
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
        total_propensity = 0;
        for (auto const &rxn : rxn_order) {
            total_propensity += rxn.propensity_val;
        }
    }

    double scNNoiSE::hill_function (int tf_count, double hill_coeff,
        double half_maximal, bool activator, double prob_contr) {
        double tf_count_pow = pow(double(tf_count), hill_coeff);
        double half_maximal_pow = pow(double(half_maximal), hill_coeff);
        if (activator) {
            return (tf_count_pow/(tf_count_pow + half_maximal_pow))*prob_contr;
        }else{
            return (half_maximal_pow/(tf_count_pow + half_maximal_pow))*prob_contr;
        }
    }

    double scNNoiSE::regulation_function (int gene_selected, int rxn) {
        double regulation_val = 0;
        for (int src = 0; src < network[0].parent_list[gene_selected].size(); ++src) {
            std::vector<int>::iterator it =
            std::find(network[0].adj_list[src].begin(),
                network[0].adj_list[src].end(),
                gene_selected);
            int gene_selected_id = std::distance(network[0].adj_list[src].begin(),
                gene_selected);
            if (network[0].edge_rxn_params[src][gene_selected_id].rxn_IN == rxn) {
                int out_species =
                    network[0].edge_rxn_params[src][gene_selected_id].species_OUT;
                double hill_coeff =
                    network[0].edge_rxn_params[src][gene_selected_id].hill_coeff;
                double half_maximal =
                    network[0].edge_rxn_params[src][gene_selected_id].half_maximal;
                double prob_contr =
                    network[0].edge_rxn_params[src][gene_selected_id].prob_contr;
                bool activator =
                    network[0].edge_rxn_params[src][gene_selected_id].activator;
                int tf_count =
                    reactions[network[0].parent_list[gene_selected][src]].molecule_count_cur[out_species];
                regulation_val +=
                    prob_contr * hill_function(tf_count, hill_coeff, half_maximal,
                                                            activator);
            }
        }
        return regulation_val;
    }


    void scNNoiSE::create_GRN_from_file (std::string filepath) {
        std::ifstream GRN_file(filepath);
        std::string row_text;
        std::vector<int> GRN_int_param;
        std::vector<double> GRN_params;
        std::vector<bool> GRN_activation;
        std::string word;
        while (std::getline(GRN_file, row_text)) {
            GRN_int_param.clear();
            GRN_params.clear();
            GRN_activation.clear();
            std::istringstream str_stream(row_text);

            unsigned int id_counter = 0;
            while (std::getline(str_stream, word, ',')) {
                if (id_counter < 2 || (id_counter >=5 && id_counter <= 6)) {
                    if (id_counter < 2) {
                        GRN_int_param.push_back(gene_rev_map[word]);
                    }else {
                        if (id_counter == 5) {
                            GRN_int_param.push_back(rxn_rev_map[word]);
                        }else {
                            GRN_int_param.push_back(species_rev_map[word]);
                        }
                    }
                }else{
                    if (id_counter == 7) {
                        std::transform(word.begin(), word.end(), word.begin(), ::tolower);
                        std::istringstream is(word);
                        bool bool_val;
                        is >> std::boolalpha >> bool_val;
                        GRN_activation.push_back(bool_val);
                    }else{
                        GRN_params.push_back(std::stod(word));
                    }
                }
                id_counter += 1;
            }
            network[0].add_edge_kinetics(GRN_int_param[0], GRN_int_param[1],
                GRN_params[0], GRN_params[1], GRN_params[2],
                GRN_int_param[2], GRN_int_param[3], GRN_activation[0]);

        }
    }
}
