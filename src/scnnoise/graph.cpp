// Graph class and algorithms
#include <vector>
#include <algorithm>
#include <iostream>
#include "graph.hpp"

namespace GRAPH {
  // Constructor
  Graph::Graph (int N) {
    // const int sz = 1000
    adj_list.reserve(N);
    parent_list.reserve(N);
    num_nodes = N;
    // adj_list.resize(N);
  }

  int Graph::get_size () {
    return num_nodes;
  }

  // // Function to add edge
  // inline void Graph::add_edge (int src, int dest, double max_expr, double hill_coeff,
  //   double half_maximal) {
  //   adj_list[src].push_back(dest);
  //   parent_list[dest].push_back(src);
  //   edge_kinetic_params[src][dest].max_expr = max_expr;
  //   edge_kinetic_params[src][dest].hill_coeff = hill_coeff;
  //   edge_kinetic_params[src][dest].half_maximal = half_maximal;
  // }

  // Function to delete edge
  void Graph::del_edge (int src, int dest) {
    adj_list[src].erase(std::remove(adj_list[src].begin(),
    adj_list[src].end(), dest), adj_list[src].end());
  }

  // Make graph undirected.
  void Graph::make_undirected (std::vector<Edge> &edges) {
    int count = 0;
    for(int i = 0; i < num_nodes; ++i){
      for (int v: adj_list[i]){
        std::vector<int>::iterator it;
        it = find(adj_list[v].begin(), adj_list[v].end(), i);
        if (it != adj_list[v].end()) {
          Edge edge;
          edge.src = v;
          edge.dest = i;
          edges.push_back(edge);
        }
        else {
          Edge edge;
          edge.src = -1;
          edge.dest = -1;
          edges.push_back(edge);
        }
        count += 1;
      }

    }

    for (auto edge: edges) {
      if (edge.src != -1) {
        adj_list[edge.src].push_back(edge.dest);
      }
    }
  }

  // Unmake graph undirected.
  void Graph::unmake_undirected (std::vector<Edge> &edges) {
    for (auto edge: edges) {
      if (edge.src != -1) {
        adj_list[edge.src].push_back(edge.dest);
        adj_list[edge.src].erase(std::remove(adj_list[edge.src].begin(),
        adj_list[edge.src].end(), edge.dest), adj_list[edge.src].end());
      }
    }
  }

  // Recursive DFS function for checking graph connectedness.
  void Graph::is_connected_DFS (int vert, std::vector<bool> visited) {
    if (!visited[vert]) {
      std::cout << "current node = " << vert << std::endl;
      visited[vert] = true;
      for (int v: adj_list[vert]) {
        if(!visited[v]){
          is_connected_DFS(v, visited);
        }
      }
    }
  }

  // Function to check whether graph is connected or not.
  bool Graph::is_connected () {
    std::vector<bool> visited(num_nodes);
    std::vector<Edge> edges;
    for (auto vert: visited) {
      vert = false;
    }
    make_undirected(edges);
    is_connected_DFS(0, visited);
    unmake_undirected(edges);
    edges.clear();
    int count = 0;
    for(int i = 0; i < num_nodes; ++i){
      if(visited[i]){
        count += 1;
      }
    }
    if(count == num_nodes){
      return true;
    }else{
      return false;
    }
  }

  // Recursive DFS utility function for determining whether graph is DAG or not.
  bool Graph::is_DAG_util (int vert, std::vector<bool> visited,
                           std::vector<bool> active) {

    visited[vert] = true;
    active[vert] = true;
    for (int v: adj_list[vert]) {
      if(!visited[v]){
        if(!is_DAG_util(v, visited, active)){
          return false;
        }
      }
      else if(active[v]){
        return false;
      }
      else{}
    }
    active[vert] = false;
    return true;
  }

  // Function to check whether graph is DAG or not.
  bool Graph::is_DAG () {
    std::vector<bool> visited(num_nodes);
    std::vector<bool> active(num_nodes);
    for (auto vert: visited) {
      vert = false;
    }
    for (auto vert: active) {
      vert = false;
    }
    for(int v = 0; v < num_nodes; ++v){
      if(!visited[v]){
        if(!is_DAG_util(v, visited, active)){
          return false;
        }
      }
    }
    return true;
  }

  // Function to return children nodes of a given node.
  void Graph::find_children (int vert, std::vector<int> &children) {
    if (!adj_list[vert].empty()) {
      for (int v: adj_list[vert]) {
        children.push_back(v);
      }
    }
  }
}
