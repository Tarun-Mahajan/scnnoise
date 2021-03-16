#include <vector>
#include <algorithm>
#include <iostream>
#include "graph.hpp"
#include "GRN.hpp"

GRN::GRN(int N):Graph(N){
}

void initialize_counts(int num_isoforms, std::vector<int> gene_active,
std::vector<int> gene_inactive, std::vector<int> nascent_mRNA,
std::vector<int> protein, std::vector<std::vector<int>> mature_mRNA){
  this.gene_active.resize(num_nodes);
  this.gene_inactive.resize(num_nodes);
  this.nascent_mRNA.resize(num_nodes);
  this.mature_mRNA.resize(num_nodes);
  this.protein.resize(num_nodes);

  this.gene_active = gene_active;
  this.gene_inactive = gene_inactive;
  this.nascent_mRNA = nascent_mRNA;
  this.mature_mRNA = mature_mRNA;
  this.protein = protein;
}
