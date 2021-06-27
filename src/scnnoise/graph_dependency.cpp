// Graph class and algorithms
#include "graph_dependency.hpp"
#include "graph.hpp"

namespace Graph_ {
  // Constructor
  GraphDependency::GraphDependency (int N): Graph_::Graph(N) {
  }

  // Function to add edge
  inline void GraphDependency::add_edge (int src, int dest) {
    adj_list[src].push_back(dest);
    parent_list[dest].push_back(src);
  }
}
