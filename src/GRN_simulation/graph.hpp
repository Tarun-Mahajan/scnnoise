// Graph header file
#ifndef GRAPH_H
#define GRAPH_H

struct Edge{
  int src;
  int dest;
};

class Graph{

  // Recursive DFS function for checking graph connectedness.
  void is_connected_DFS(int vert, bool visited[]);

  // Make graph undirected.
  void make_undirected(std::vector<Edge> &edges);

  // Unmake graph undirected.
  void unmake_undirected(std::vector<Edge> &edges);

  // Recursive DFS utility function for determining whether graph is DAG or not.
  bool is_DAG_util(int vert, bool visited[], bool active[]);



public:
  // Number of nodes in the Graph.
  int num_nodes;

  // A vector or vectors to represent adjacency list.
  std::vector<std::vector<int>> adj_list;

  // Graph constructor.
  Graph(int N);

  // Function to add edge.
  void add_edge(int src, int dest);

  // Function to delete edge.
  void del_edge(int src, int dest);

  // Function to check whether graph is connected or not.
  bool is_connected();

  // Function to check whether graph is DAG or not
  bool is_DAG();
};

#endif
