// GRN header
#ifndef GRN_H
#define GRN_H

class GRN: public Graph{
  std::vector<int> gene_active, gene_inactive, nascent_mRNA, protein;
  std::vector<std::vector<int>> mature_mRNA;

public:
  GRN(int N);
  void initialize_counts(int num_isoforms);

};


#endif
