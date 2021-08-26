from scnnoise import _scnnoise
import pandas as pd

class LTSim:
    def __init__(self, root, num_rxns, num_genes, num_species_genes, num_rxns_gene_type, max_time, count_save_file, sample_csv, edge_list_csv, sample_species):
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
        self.sample_csv= sample_csv
        self.edge_list_csv = edge_list_csv
        self.sample_species = sample_species
        self.count_save_file = count_save_file
        
    def sim_LT(self, num_samples):
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
        
        gene_type = [0,0] 
        molecule_count_cur = [[0,2,0,0],[0,2,0,0]] #initialize based on self.rxn_rate

        GRN_rxn_IN = [[], [2]] #Initialize w/ csv?
        GRN_species_OUT = [[3], []] #Initialize w/ csv?
        

        #This should all be based on gene type str
        reactants = [[[1],[0],[0],[2],[2],[3]], [[1],[0],[0],[2],[2],[3]]]
        products = [[[0],[1],[0,2],[],[2,3],[]], [[0],[1],[0,2],[],[2,3],[]]]
        reactants_stoichio = [[[1],[1],[1],[1],[1],[1]], [[1],[1],[1],[1],[1],[1]]]
        products_stoichio = [[[1],[1],[1,1],[],[1,1],[]], [[1],[1],[1,1],[],[1,1],[]]]
        
        propensity_val = [[2,0,0,0,0,0], [2,0,0,0,0,0]] #calculated internally?


        for gene_id in range(self.num_genes):
            simulator.add_gene_state(gene_id, gene_type[gene_id] ,GRN_rxn_IN[gene_id], GRN_species_OUT[gene_id], molecule_count_cur[gene_id], reactants[gene_id], products[gene_id], reactants_stoichio[gene_id], products_stoichio[gene_id], self.root.rxn_rate[gene_id], propensity_val[gene_id])
        
    #void scNNoiSE::add_GRN_edge (int src, int dest, double prob_contr,
    #double hill_coeff, double half_maximal, int rxn_IN, int species_OUT,
    #bool activator)
        
        # Use file name for this
        edge_list = pd.read_csv(self.edge_list_csv)
        for edge in edge_list.values:
            simulator.add_GRN_edge(edge[0], edge[1], edge[2], edge[3], edge[4], edge[5],edge[6], edge[7])
    
        #this should be set up based on gene types in gene_type
        simulator.add_dependency_edge(0,0,1)
        simulator.add_dependency_edge(0,0,2)
        simulator.add_dependency_edge(0,1,0)
        simulator.add_dependency_edge(0,1,2)
        simulator.add_dependency_edge(0,2,3)
        simulator.add_dependency_edge(0,3,4)
        simulator.add_dependency_edge(0,4,5)
    
        
        #Step 2: Run sim_cell_type() on root
        self.root.sim_transition(num_samples, simulator, False)
        
        
        #Step 3: return result file, or data structure
            
            
        return pd.read_csv(self.sample_csv)