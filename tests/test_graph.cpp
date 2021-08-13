#include <iostream>
#include <algorithm>
#include "../src/scnnoise/graph_derived.hpp"

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

int main () {
    // Set graph edge parameters
    int num_nodes = 6;
    double prob_contr = 0.1;
    double hill_coeff = 2;
    double half_maximal = 10;
    int rxn_IN = 0;
    int species_OUT = 3;
    bool activator = true;

    // Create graph
    GraphSpace::GRN gene_net(num_nodes);
    gene_net.add_edge_kinetics(0, 1, prob_contr, hill_coeff, half_maximal,
                            rxn_IN, species_OUT, activator);
    gene_net.add_edge_kinetics(1, 2, prob_contr, hill_coeff, half_maximal,
                            rxn_IN, species_OUT, activator);
    gene_net.add_edge_kinetics(2, 3, prob_contr, hill_coeff, half_maximal,
                            rxn_IN, species_OUT, activator);
    gene_net.add_edge_kinetics(3, 4, prob_contr, hill_coeff, half_maximal,
                            rxn_IN, species_OUT, activator);
    gene_net.add_edge_kinetics(0, 4, prob_contr, hill_coeff, half_maximal,
                            rxn_IN, species_OUT, activator);
    gene_net.add_edge_kinetics(4, 0, prob_contr, 2*hill_coeff, half_maximal,
                            rxn_IN, species_OUT, activator);
    gene_net.add_edge_kinetics(1, 5, prob_contr, 3*hill_coeff, half_maximal,
                            rxn_IN, species_OUT, activator);

    // Create testing flags
    bool is_num_nodes_set = false;
    bool is_count_edges_correct = false;
    bool is_adjlist_set = false;
    bool is_parentlist_set = false;
    bool are_edgeparams_set = false;
    bool is_find_children_working = false;
    bool is_find_children_edge_info_working = false;
    unsigned int num_edges = 7;

    if (count_edges(gene_net) == num_edges) {
        is_count_edges_correct = true;
    }

    if (gene_net.num_nodes == num_nodes) {
        is_num_nodes_set = true;
    }

    // Test adj_list
    unsigned int count_adj = 0;
    // check 0-->1
    int src = 0;
    int dest = 1;
    std::vector<int>::iterator it = std::find(gene_net.adj_list[src].begin(),
        gene_net.adj_list[src].end(), dest);
    if (it != gene_net.adj_list[src].end()) {
        ++count_adj;
    }

    // check 1-->2
    src = 1;
    dest = 2;
    it = std::find(gene_net.adj_list[src].begin(),
        gene_net.adj_list[src].end(), dest);
    if (it != gene_net.adj_list[src].end()) {
        ++count_adj;
    }

    // check 2-->3
    src = 2;
    dest = 3;
    it = std::find(gene_net.adj_list[src].begin(),
        gene_net.adj_list[src].end(), dest);
    if (it != gene_net.adj_list[src].end()) {
        ++count_adj;
    }

    // check 3-->4
    src = 3;
    dest = 4;
    it = std::find(gene_net.adj_list[src].begin(),
        gene_net.adj_list[src].end(), dest);
    if (it != gene_net.adj_list[src].end()) {
        ++count_adj;
    }

    // check 0-->4 and 4-->0
    src = 0;
    dest = 4;
    it = std::find(gene_net.adj_list[src].begin(),
        gene_net.adj_list[src].end(), dest);
    if (it != gene_net.adj_list[src].end()) {
        ++count_adj;
    }

    src = 4;
    dest = 0;
    it = std::find(gene_net.adj_list[src].begin(),
        gene_net.adj_list[src].end(), dest);
    if (it != gene_net.adj_list[src].end()) {
        ++count_adj;
    }

    src = 1;
    dest = 5;
    it = std::find(gene_net.adj_list[src].begin(),
        gene_net.adj_list[src].end(), dest);
    if (it != gene_net.adj_list[src].end()) {
        ++count_adj;
    }

    if (count_adj == num_edges) {
        is_adjlist_set = true;
    }

    // Test parent_list
    unsigned int count_parent = 0;
    // check 0-->1
    src = 0;
    dest = 1;
    it = std::find(gene_net.parent_list[dest].begin(),
        gene_net.parent_list[dest].end(), src);
    if (it != gene_net.parent_list[dest].end()) {
        ++count_parent;
    }

    // check 1-->2
    src = 1;
    dest = 2;
    it = std::find(gene_net.parent_list[dest].begin(),
        gene_net.parent_list[dest].end(), src);
    if (it != gene_net.parent_list[dest].end()) {
        ++count_parent;
    }

    // check 2-->3
    src = 2;
    dest = 3;
    it = std::find(gene_net.parent_list[dest].begin(),
        gene_net.parent_list[dest].end(), src);
    if (it != gene_net.parent_list[dest].end()) {
        ++count_parent;
    }

    // check 3-->4
    src = 3;
    dest = 4;
    it = std::find(gene_net.parent_list[dest].begin(),
        gene_net.parent_list[dest].end(), src);
    if (it != gene_net.parent_list[dest].end()) {
        ++count_parent;
    }

    // check 0-->4 and 4-->0
    src = 0;
    dest = 4;
    it = std::find(gene_net.parent_list[dest].begin(),
        gene_net.parent_list[dest].end(), src);
    if (it != gene_net.parent_list[dest].end()) {
        ++count_parent;
    }

    src = 4;
    dest = 0;
    it = std::find(gene_net.parent_list[dest].begin(),
        gene_net.parent_list[dest].end(), src);
    if (it != gene_net.parent_list[dest].end()) {
        ++count_parent;
    }

    src = 1;
    dest = 5;
    it = std::find(gene_net.parent_list[dest].begin(),
        gene_net.parent_list[dest].end(), src);
    if (it != gene_net.parent_list[dest].end()) {
        ++count_parent;
    }

    if (count_adj == num_edges) {
        is_parentlist_set = true;
    }

    // Test find children
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


    std::cout << "Num of edges correct = " << is_count_edges_correct << std::endl;
    std::cout << "Num of nodes properly set = " << is_num_nodes_set << std::endl;
    std::cout << "Adjacency list properly set = " << is_adjlist_set << std::endl;
    std::cout << "Parent list properly set = " << is_parentlist_set << std::endl;
    std::cout << "find_children properly working = " << is_find_children_working << std::endl;

    return 0;
}
