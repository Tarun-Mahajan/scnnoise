from scnnoise import _scnnoise
import pandas as pd
import numpy as np

class CellType:
    def __init__(self, lineageName,rxn_rates, children):
        """
        Takes in rxn_rates and children and constructs a cellType node
        """
        #get rate parameters from external file or otherwise
        self.lineageName = lineageName
        self.rxn_rates = rxn_rates
        self.children = children
        self.count_csv = count_csv
        self.sample_csv = sample_csv

    def sim_cell_type(num_samples, simulator):
        """
        Params
            num_samples - number of sample reads to output for this cell type
            simulator -  a simulator object used to simulate  the system
            
        Simulates num_samples reads of cellType 
        and outputs num_samples x genes Matrix of raw mRNA counts
        """
        #Recursive function (careful with passing by reference/by copy)
        #Step 1: simulate gene expression at steady state
            #a. it would be nice to be able to pass the 4 to simulate for into the simulate function
        f = open(self.count_csv, "w")
        f.truncate()
        f.close()
        simulator.simulate(1000)

            #for steady state check sergio method
        
        #Step 2: Sample num_samples reads for each gene
        sim_out = pd.read_csv(self.count_csv)
        samples = np.random.randint(0,len(sim_out.index), size = num_samples)
        sample_out = sim_out.iloc[samples]
        sample_out['Cell Type'] = [str(self.LineageName)] * num_samples
        sample_out.to_csv(self.sample_csv, mode = 'a')

        #Step 3: Recusively run sim_transition() on all children (this could be parallelized)
        #recursive case
        if len(children) != 0:
            for cell_type in self.children:
                cell_type.sim_transition(num_samples, simulator)
                
        
    def sim_transition(num_samples, simulator):
        """
        Params
            num_samples - number of sample reads to output for this cell type
            simulator -  a simulator object used to simulate  the system
            
        Simulates num_samples reads of cells transitioning from the parent cell type 
        and outputs num_samples x genes Matrix of raw mRNA counts
        """
        #helper (careful with passing by reference/by copy)
        #Step 1: Save new kinetic parameters into simulator
            #a. maybe here it would be better to add options to the c++ code to make this easier 
        for gene in range(simulator.num_genes):   
            simulator.set_rxn_rates(self.rxn_rates)
            
            
        #Step 2: simluate transition from parent steady state to daughter steady state
            #a. here we check when the system reaches a new steady state and sample before then
            #b. if steady states are similar mbe set simulation time?
            
        #try t-test for mean for all genes (for each thousand timepoints) compared to calculated steady state for current cell_type
        f = open(self.count_csv, "w")
        f.truncate()
        f.close()        
        simulator.simulate(1000)
        sim_out = pd.read_csv(self.count_csv)
        #identify molecules of interest in csv

        #steady state detection

        #collect Traisition samples
        samples = np.random.randint(0,len(sim_out.index), size = num_samples)
        sample_out = sim_out.iloc[samples]
        sample_out['Cell Type'] = [str(self.LineageName)+'T'] * num_samples
        sample_out.to_csv(self.sample_csv, mode = 'a')
        
        #Step 3: Sample num_sample transition cells and store to output and run sim_cell_type
        self.sim_cell_type(num_samples, simulator)
        