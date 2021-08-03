from scnnoise import _scnnoise
import pandas as pd

class LTSim:
    def __init__(self, root, num_rxns, num_genes, num_species_genes, num_rxns_gene_type, max_time, count_save_file, edge_list_csv):
        """
        Constructor for the LTSim object
        Takes in the root of the LineageTree in the simulation
        """
        self.root = root    #cell Type
        
        #write set functions for each instead of in constructor
        self.num_rxns = num_rxns
        self.num_genes = num_genes
        self.num_species_genes = num_species_genes
        self.num_rxns_gene_type = num_rxns_gene_type
        self.max_time = max_time
        self.count_save_file = count_save_file
        self.edge_list_csv = edge_list_csv
        
    def sim_LT(num_samples):
        """
        Outputs a n x g Matrix (n cells by g Genes) of raw mRNA counts
        """
        #Step 1: Construct simulator object
        #a. build GRN from params
        #b. define intragene interactions
        simulator = _scnnoise.gillespieSDMnoCellCycle(self.num_rxns, self.num_genes, self.num_species_genes, self.num_rxns_gene_type, self.max_time, True, 10000, self.count_save_file)
        
    #add_gene_state (int gene_id, int gene_type, std::vector<int> GRN_rxn_IN,
    #std::vector<int> GRN_species_OUT, std::vector<int> molecule_count_cur,
    #std::vector<std::vector<int>> reactants, std::vector<std::vector<int>> products,
    #std::vector<std::vector<int>> reactants_stoichio, std::vector<std::vector<int>> products_stoichio,
    #std::vector<double> rxn_rate, std::vector<double> propensity_val)
        
        for gene_id in range(self.num_genes):
            simulator.add_gene_state(gene_id, gene_type[gene_id] ,GRN_rxn_IN[gene_id], GRN_species_OUT[gene_id], molecule_count_cur[gene_id], reactants[gene_id], products[gene_id], reactants_stoichio[gene_id], products_stoichio[gene_id], rxn_rate[gene_id], products_stoichio[gene_id], rxn_rate[gene_id], propensity_val[gene_id])
        
    #void scNNoiSE::add_GRN_edge (int src, int dest, double prob_contr,
    #double hill_coeff, double half_maximal, int rxn_IN, int species_OUT,
    #bool activator)
        
        # Use file name for this
        for edge in edge_list:
            simulator.add_GRN_edge(edge[0], edge[1], edge[2], edge[3], edge[4], edge[5],edge[6], edge[7])
        
        
        for dep in dep_list:
            simulator.add_dependency_edge(dep[0], dep[1], dep[2])

        
        #Step 2: Run sim_cell_type() on root
        self.root.sim_cell_type(num_samples, simulator)
        
        
        #Step 3: return result file, or data structure
            
            
        return self.sample