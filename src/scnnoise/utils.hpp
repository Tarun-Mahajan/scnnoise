#ifndef UTILS_H
#define UTILS_H

#include "../src/scnnoise/graph_derived.hpp"

bool test_node_count (GraphSpace::GRN gene_net, unsigned int num_nodes);

void create_GRN_from_file (GraphSpace::GRN &gene_net, std::string filepath);

unsigned int count_edges (GraphSpace::GRN gene_net);

bool test_edges_count (GraphSpace::GRN gene_net, unsigned int num_nodes);

unsigned int count_GRN_edges_from_file (std::string filepath);

bool test_find_children (GraphSpace::GRN gene_net);

std::vector<bool> test_edge_properties (GraphSpace::GRN gene_net, std::string filepath)


#endif
