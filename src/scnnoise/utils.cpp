#include <iostream>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <string>
#include "../src/scnnoise/graph_derived.hpp"
#include "utils.hpp"

void create_GRN_from_file (GraphSpace::GRN &gene_net, std::string filepath,
    unsigned int num_nodes) {
    std::ifstream GRN_file(filepath);
    std::string row_text;
    vector<int> GRN_int_param;
    vector<double> GRN_params;
    vector<bool> GRN_activation;
    std::string word;
    while (std::getline(GRN_file, row_text)) {
        GRN_int_param.clear();
        GRN_params.clear();
        GRN_activation.clear();
        std::stringstream str_stream(row_text);

        unsigned int id_counter = 0;
        while (std::getline(str_stream, word, ', ')) {
            if (id_counter < 4) {
                GRN_int_param.push_back(std::stoi(word));
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
        }
        gene_net.add_edge_kinetics(GRN_int_param[0], GRN_int_param[1],
            GRN_params[0], GRN_params[1], GRN_params[2],
            GRN_int_param[2], GRN_int_param[3], GRN_activation[0]);

    }
}

unsigned int count_edges (GraphSpace::GRN gene_net) {
    unsigned int num_edges = 0;
    unsigned int node = 0;
    for (auto src = gene_net.adj_list.begin(); src != gene_net.adj_list.end(); ++src) {
        for (auto dest = (*src).begin(); dest != (*src).end(); ++ dest) {
        }
        num_edges += (*src).size();
        node += 1;
    }
    return num_edges;
}

unsigned int count_GRN_edges_from_file (std::string filepath) {
    std::ifstream GRN_file(filepath);
    std::string row_text;
    unsigned int count_edges = 0;
    while (std::getline(GRN_file, row_text)) {
        count_edges += 1;
    }
    return count_edges;
}

bool test_edges_count (GraphSpace::GRN gene_net, std::string filepath) {
    unsigned int count_edges_from_file =
        count_GRN_edges_from_file(std::string filepath);
    unsigned int count_edges_GRN = count_edges(gene_net);
    if (count_edges_from_file == count_edges_GRN) {
        return true;
    }else{
        return false;
    }
}

bool test_node_count (GraphSpace::GRN gene_net, unsigned int num_nodes) {
    if (gene_net.num_nodes == num_nodes) {
        return true;
    }else{
        return false;
    }
}


std::vector<bool> test_edge_properties (GraphSpace::GRN gene_net, std::string filepath) {
    std::vector<bool> out_flags;
    bool is_adjlist_set = true;
    bool is_parentlist_set = true;
    bool are_edgeparams_set = false;
    std::ifstream GRN_file(filepath);
    std::string row_text;
    vector<int> GRN_int_param;
    vector<double> GRN_params;
    vector<bool> GRN_activation;
    std::string word;
    while (!(is_adjlist_set == false && is_parentlist_set == false &&
        are_edgeparams_set == false)) {
        while (std::getline(GRN_file, row_text)) {
            GRN_int_param.clear();
            GRN_params.clear();
            GRN_activation.clear();
            std::stringstream str_stream(row_text);

            unsigned int id_counter = 0;
            while (std::getline(str_stream, word, ', ')) {
                if (id_counter < 4) {
                    GRN_int_param.push_back(std::stoi(word));
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
            }
            std::vector<int>::iterator it =
                std::find(gene_net.adj_list[GRN_int_param[0]].begin(),
                    gene_net.adj_list[GRN_int_param[0]].end(), GRN_int_param[1]);
            if (is_adjlist_set) {
                if (it == gene_net.adj_list[GRN_int_param[0]].end()) {
                    is_adjlist_set = false;
                    are_edgeparams_set = false;
                }else{
                    if (are_edgeparams_set) {
                        int id =
                            std::distance(gene_net.adj_list[GRN_int_param[0]].begin(), GRN_int_param[1]);
                        unsigned int count_param_test = 0;
                        if (gene_net.edge_rxn_params[GRN_int_param[0]][id].prob_contr ==
                            GRN_params[0]) {
                                count_param_test += 1;
                            }

                        if (gene_net.edge_rxn_params[GRN_int_param[0]][id].hill_coeff ==
                            GRN_params[1]) {
                                count_param_test += 1;
                            }

                        if (gene_net.edge_rxn_params[GRN_int_param[0]][id].half_maximal ==
                            GRN_params[2]) {
                                count_param_test += 1;
                            }

                        if (gene_net.edge_rxn_params[GRN_int_param[0]][id].rxn_IN ==
                            GRN_int_param[2]) {
                                count_param_test += 1;
                            }

                        if (gene_net.edge_rxn_params[GRN_int_param[0]][id].species_OUT ==
                            GRN_int_param[3]) {
                                count_param_test += 1;
                            }

                        if (gene_net.edge_rxn_params[GRN_int_param[0]][id].activator ==
                            GRN_activation[0]) {
                                count_param_test += 1;
                            }

                        if (count_param_test != 6) {
                            are_edgeparams_set = false;
                        }
                    }
                }
            }

            it = std::find(gene_net.parent_list[GRN_int_param[1]].begin(),
                    gene_net.parent_list[GRN_int_param[1]].end(), GRN_int_param[0]);
            if (is_parentlist_set) {
                if (it == gene_net.parent_list[GRN_int_param[1]].end()) {
                    is_parentlist_set = false;
                }
            }

        }
    }
    out_flags.push_back(is_adjlist_set);
    out_flags.push_back(is_parentlist_set);
    out_flags.push_back(are_edgeparams_set);
    return out_flags;
}

bool test_find_children (GraphSpace::GRN gene_net) {
    bool is_find_children_working = false;
    unsigned int num_nodes = gene_net.num_nodes;
    std::vector<int> children;
    unsigned int num_correct_children_nodes = 0;

    for (unsigned int src = 0; src < num_nodes; ++src) {
        gene_net.find_children (src, children);
        unsigned int num_correct_children = 0;
        for (auto dest = gene_net.adj_list[src].begin();
            dest != gene_net.adj_list[src].end(); ++dest) {
                std::vector<int>::iterator it = std::find(children.begin(),
                    children.end(), *dest);
                if (it != children.end()) {
                    num_correct_children += 1;
                }
            }
        if ((num_correct_children == gene_net.adj_list[src].size()) &&
            (children.size() == gene_net.adj_list[src].size())) {
                num_correct_children_nodes += 1;
            }
    }
    if (num_correct_children_nodes == num_nodes) {
        is_find_children_working = true;
    }
    return is_find_children_working;
}
