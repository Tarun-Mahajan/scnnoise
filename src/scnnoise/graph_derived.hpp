// Graph header file
#ifndef GRN_H
#define GRN_H

#include "graph.hpp"
#include <vector>

struct edge_rxn_struct {
  double prob_contr;
  double hill_coeff;
  double half_maximal;
  int rxn_IN;
  int species_OUT;
  bool activator;
};

namespace GraphSpace {
  class GraphDependency : public Graph  {

  public:
    /* Memeber functions */
    /********************************************//**
     \brief Graph Constructor.

     \param[in] N Number of nodes in the graph.
     ***********************************************/
    GraphDependency (int N);


    /********************************************//**
     \brief Function to add edge to the graph.

     \param[in] src source vertex for the edge.
     \param[in] dest destination vertex for the edge.
     ***********************************************/
    void add_edge (int src, int dest) override;
  };

  // class GRN
  class GRN : public Graph  {
  protected:
    std::vector<std::vector<edge_rxn_struct>> edge_rxn_params;

  public:
    /* Memeber functions */
    /********************************************//**
     \brief Graph Constructor.

     \param[in] N Number of nodes in the graph.
     ***********************************************/
    GRN (int N);



    void add_edge_kinetics (int src, int dest, double prob_contr,
                            double hill_coeff, double half_maximal,
                            int rxn_IN, int species_OUT);
    /********************************************//**
     \brief Function to add edge to the graph.

     \param[in] src source vertex for the edge.
     \param[in] dest destination vertex for the edge.
     ***********************************************/
    void add_edge (int src, int dest) override;

    void find_children_edge_info (int vert, std::vector<edge_rxn_struct> &children_edge_info);
  };
}

#endif
