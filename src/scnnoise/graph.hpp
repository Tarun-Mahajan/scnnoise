// Graph header file
#ifndef GRAPH_H
#define GRAPH_H

#include <vector>

struct Edge {
  int src;
  int dest;
};

namespace GraphSpace {
    class Graph {
    public:
        /* Data */
        /********************************************//**
        \brief Number of nodes in the graph.
        ***********************************************/
        int num_nodes;

        /********************************************//**
        \brief A vector of vectors to represent adjacency list.
        ***********************************************/
        std::vector<std::vector<int>> adj_list;
        std::vector<std::vector<int>> parent_list;
        // std::vector<std::map<int, edge_struct>> edge_kinetic_params;

        /* Memeber functions */
        /********************************************//**
        \brief Recursive DFS function for checking graph connectedness.

        \param[in] vert vertex in the graph from where the DFS is initiated.
        \param[in] visited boolean vector to track whether vertices in the
                graph have been visited by the DFS or not.
        ***********************************************/
        void is_connected_DFS (int vert, std::vector<bool> visited);

        // Make graph undirected.
        /********************************************//**
        \brief Make graph undirected for checking weak connectedness.

        \param[in] vector of edges to store edges added to make the graph
                undirected.
        ***********************************************/
        void make_undirected (std::vector<Edge> &edges);

        /********************************************//**
        \brief Unmake graph undirected after checking for weak connectedness.

        \param[in] vector of edges to store edges to delete
                to make the graph directed again.
        ***********************************************/
        void unmake_undirected (std::vector<Edge> &edges);

        // Recursive DFS utility function for determining whether graph is DAG or not.
        /********************************************//**
        \brief Recursive DFS utility function for determining whether
            graph is DAG or not.

        \param[in] vert vertex in the graph from where the DFS is initiated.
        \param[in] visited boolean vector to track whether vertices in the
                graph have been visited by the DFS or not.
        ***********************************************/
        bool is_DAG_util (int vert, std::vector<bool> visited,
            std::vector<bool> active);



        // public:
        /* Memeber functions */
        /********************************************//**
        \brief Graph Constructor.

        \param[in] N Number of nodes in the graph.
        ***********************************************/
        Graph (int N);

        int get_size ();

        /********************************************//**
        \brief Function to add edge to the graph.

        \param[in] src source vertex for the edge.
        \param[in] dest destination vertex for the edge.
        ***********************************************/
        virtual void add_edge (int src, int dest) = 0;

        /********************************************//**
        \brief Function to delete edge from the graph.

        \param[in] src source vertex for the edge.
        \param[in] dest destination vertex for the edge.
        ***********************************************/
        void del_edge (int src, int dest);

        /********************************************//**
        \brief Function to check whether graph is weakly connected or not.
        ***********************************************/
        bool is_connected ();

        /********************************************//**
        \brief Function to check whether graph is DAG or not.
        ***********************************************/
        bool is_DAG ();

        /********************************************//**
        \brief Function to return children nodes of a given node.

        \param[in] vert vertex to return children for.
        \param[in] children vector to store the children nodes.
        ***********************************************/
        void find_children (int vert, std::vector<int> &children);

        void find_parents (int vert, std::vector<int> &parents);
    };
}

#endif
