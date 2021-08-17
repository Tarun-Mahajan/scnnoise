#ifndef UTILS_H
#define UTILS_H

#include <string>
#include "graph_derived.hpp"

void create_GRN_from_file (GraphSpace::GRN &gene_net, std::string filepath,
    std::map<std::string, int> gene_rev_map);

#endif
