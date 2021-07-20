class LTSim:
    def __init__(self, root):
        """
        Constructor for the LTSim object
        Takes in the root of the LineageTree in the simulation
        """
        self.root = root    #cell Type

    def sim_LT(num_samples):
        """
        Outputs a n x g Matrix (n cells by g Genes) of raw mRNA counts
        """
        #Step 1: Construct simulator object
        #a. build GRN from params
        #b. define intragene interactions
        
        #Step 2: Run sim_cell_type() on root
        self.root.sim_cell_type(num_samples, simulator)
        
        #Step 3: return result file, or data structure
            
        return self.sample