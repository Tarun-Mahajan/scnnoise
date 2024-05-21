/* Copyright (C) 2021 Tarun Mahajan - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GNU AGPLv3 license.
 *
 * You should have received a copy of the GNU AGPLv3 license with
 * this file. If not, please write to: tarunm3@illinois.edu, or 
 * visit : https://www.gnu.org/licenses/agpl-3.0.en.html
 */

#ifndef UTILS_H
#define UTILS_H

#include <string>
#include "graph_derived.hpp"

void create_GRN_from_file (GraphSpace::GRN &gene_net, std::string filepath,
    std::map<std::string, int> gene_rev_map);

#endif
