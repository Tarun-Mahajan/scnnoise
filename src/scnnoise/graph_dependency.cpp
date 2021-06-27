// Graph class and algorithms
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include "graph.hpp"

// Constructor
GraphDependency::GraphDependency (int N): Graph(N) {
}

// Function to add edge
inline void GraphDependency::add_edge (int src, int dest) {
  adj_list[src].push_back(dest);
  parent_list[dest].push_back(src);
}
