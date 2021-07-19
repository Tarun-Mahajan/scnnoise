class cellType:
    def __init__(self, lineageName,rxn_rates, children):
        """
        Takes in rxn_rates and children and constructs a cellType node
        """
        #get rate parameters from external file or otherwise
        self.lineageName = lineageName
        self.rxn_rates = rxn_rates
        self.children = children

    def sim_cell_type(num_samples, simulator):
        """
        Params
            num_samples - number of sample reads to output for this cell type
            simulator -  a simulator object used to simulate  the system
            
        Simulates num_samples reads of cellType 
        and outputs num_samples x genes Matrix of raw mRNA counts
        """
        #Recursive function (careful with passing by reference/by copy)
        #Step 1: simulate gene expression to steady state
        #Step 2: Sample num_samples reads for each gene
        #Step 3: Recusively run sim_transition() and sim_cell_type() on all children (this could be parallelized)
        
        #base case
        if len(children) == 0:
            return 
        
        #recursive case
        else:
            for cell_type in self.children:
                cell_type.sim_transition(num_samples, simulator)
                cell_type.sim_cell_type(num_samples, simulator)
        
    def sim_transition(num_samples, simulator):
        """
        Params
            num_samples - number of sample reads to output for this cell type
            simulator -  a simulator object used to simulate  the system
            
        Simulates num_samples reads of cells transitioning from the parent cell type 
        and outputs num_samples x genes Matrix of raw mRNA counts
        """
        #helper (careful with passing by reference/by copy)
        #Step 1: simluate transition from parent steady state to daughter steady state
            #a. here we check when the system reaches a new steady state and sample before then
            #b. if steady states are similar mbe set simulation time?
        #Step 2: Sample num_sample transition cells and store to output
        
        