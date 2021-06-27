// Graph header file
#ifndef GRN_H
#define GRN_H

#include <vector>

struct edge_rxn_struct {
  double prob_contr;
  double hill_coeff;
  double half_maximal;
  int rxn_IN;
  int species_OUT;
};

namespace Graph_ {
  namespace GraphDependency_ {
    class GraphDependency: public Graph  {

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
      void add_edge (int src, int dest);
    };
  }
}

#endif
