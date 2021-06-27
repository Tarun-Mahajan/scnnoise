// Graph header file
#ifndef GRN_H
#define GRN_H

#include <vector>
#include "graph.hpp"

struct edge_rxn_struct {
  double prob_contr;
  double hill_coeff;
  double half_maximal;
  int rxn_IN;
  int species_OUT;
};

namespace Graph_ {
  namespace GraphDerived {
    class GRN: public Graph  {
    protected:
      std::vector<std::vector<edge_rxn_struct>> edge_rxn_params;

    public:
      /* Memeber functions */
      /********************************************//**
       \brief Graph Constructor.

       \param[in] N Number of nodes in the graph.
       ***********************************************/
      GRN (int N);


      /********************************************//**
       \brief Function to add edge to the graph.

       \param[in] src source vertex for the edge.
       \param[in] dest destination vertex for the edge.
       ***********************************************/
      void add_edge (int src, int dest, double prob_contr, double hill_coeff,
        double half_maximal, int rxn_IN, int species_OUT);

      void find_children_edge_info (int vert, std::vector<edge_rxn_struct> &children_edge_info);
    };
  }
}

#endif
