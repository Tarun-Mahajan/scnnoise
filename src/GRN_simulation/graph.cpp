// Graph class and algorithms
#include <vector>
#include <algorithm>
#include <iostream>
#include "graph.hpp"

Graph::Graph(int N){
  num_nodes = N;
  adj_list.resize(N);
}

void Graph::add_edge(int src, int dest){
  adj_list[src].push_back(dest);
}

void Graph::del_edge(int src, int dest){
  adj_list[src].erase(std::remove(adj_list[src].begin(),
  adj_list[src].end(), dest), adj_list[src].end());
}

void Graph::make_undirected(std::vector<Edge> &edges){
  int count = 0;
  for(int i = 0; i < num_nodes; ++i){
    for (int v: adj_list[i]){
      std::vector<int>::iterator it;
      it = find(adj_list[v].begin(), adj_list[v].end(), i);
      if(it != adj_list[v].end()){
        Edge edge;
        edge.src = v;
        edge.dest = i;
        edges.push_back(edge);
      }
      else{
        Edge edge;
        edge.src = -1;
        edge.dest = -1;
        edges.push_back(edge);
      }
      count += 1;
    }

  }

  for(auto &edge: edges){
    if(edge.src != -1){
      adj_list[edge.src].push_back(edge.dest);
    }
  }
}

void Graph::unmake_undirected(std::vector<Edge> &edges){
  for(auto &edge: edges){
    if(edge.src != -1){
      adj_list[edge.src].push_back(edge.dest);
      adj_list[edge.src].erase(std::remove(adj_list[edge.src].begin(),
      adj_list[edge.src].end(), edge.dest), adj_list[edge.src].end());
    }
  }
}

void Graph::is_connected_DFS(int vert, bool visited[]){
  if(!visited[vert]){
    std::cout << "current node = " << vert << std::endl;
    visited[vert] = true;
    for (int v: adj_list[vert]){
      if(!visited[v]){
        is_connected_DFS(v, visited);
      }
    }
  }
}

bool Graph::is_connected(){
  bool* visited = new bool[num_nodes];
  std::vector<Edge> edges;
  for(int i = 0; i < num_nodes; ++i ){
    visited[i] = false;
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
  delete [] visited;
  if(count == num_nodes){
    return true;
  }else{
    return false;
  }
}

bool Graph::is_DAG_util(int vert, bool visited[], bool active[]){

  visited[vert] = true;
  active[vert] = true;
  for (int v: adj_list[vert]){
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

bool Graph::is_DAG(){
  bool* visited = new bool[num_nodes];
  bool* active = new bool[num_nodes];
  for(int i = 0; i < num_nodes; ++i){
    visited[i] = false;
    active[i] = false;
  }

  for(int v = 0; v < num_nodes; ++v){
    if(!visited[v]){
      if(!is_DAG_util(v, visited, active)){
        delete [] visited;
        delete [] active;
        return false;
      }
    }
  }
  delete [] visited;
  delete [] active;
  return true;
}
