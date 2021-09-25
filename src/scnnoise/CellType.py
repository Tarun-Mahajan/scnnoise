#SIMULATOR NEEDS DEEPCOPY BADLY
from scnnoise import _scnnoise
import pandas as pd
import numpy as np

class CellType:
    def __init__(self, lineageName,rxn_rate_csv, children):
        """
        Takes in rxn_rates and children and constructs a cellType node
        """
        #get rate parameters from external file or otherwise
        self.lineageName = lineageName
        self.children = children
        cell_types = pd.read_csv(rxn_rate_csv)
        curr_cell_type = cell_types[cell_types['Cell Type'] == self.lineageName]
        curr_cell_type.index = curr_cell_type['Gene']
        curr_cell_type = curr_cell_type.drop(columns = ['Cell Type', 'Gene'])
        self.rxn_rates = curr_cell_type.T.to_dict()

    def sim_cell_type(self, num_samples, simulator, molecule, count_csv, sample_csv):
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
        print(str(self.lineageName))
        f = open(count_csv, "w")
        f.truncate()
        f.close()
        simulator.set_simulation_params(1000, True)
        simulator.simulate(list(np.random.randint(low = 1, high = 1000000, size = 4)))

            #for steady state check sergio method
        
        #Step 2: Sample num_samples reads for each gene
        sim_out = pd.read_csv(count_csv)
        sim_out = sim_out.loc[:, ~sim_out.columns.str.contains('^Unnamed')]
        species = [col for col in sim_out.columns if col[-len(molecule):] == molecule]
        samples = np.random.random(size = num_samples) *1000
        time = np.cumsum(sim_out[sim_out.columns[0]])
        time.append(np.inf)
        idxs = np.searchsorted(time, samples,side = 'right') - 1
        sample_out = sim_out.iloc[idxs][species]
        sample_out['Cell Type'] = [str(self.lineageName)] * num_samples
        sample_out.to_csv(sample_csv, mode = 'a', header =False)

        genes = [gene.split(':') for gene in sim_out.columns[1:]]
        init_mol_count = { gene[0] : {} for gene in genes}

        for gene in genes:
            init_mol_count[ gene[0] ].update( {gene[1]: int(sim_out[ gene[0]+':'+gene[1] ].iloc[-1]) })


        #Step 3: Recusively run sim_transition() on all children (this could be parallelized)
        #recursive case
        if len(self.children) != 0:
            for cell_type in self.children:
                cell_type.sim_transition(num_samples, simulator, True, molecule, count_csv, sample_csv, init_mol_count)
                
        
    def sim_transition(self, num_samples, simulator, collect_samples, molecule, count_csv, sample_csv, init_mol_count  = {}):
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
        print(str(self.lineageName)+'T')
        simulator.swap_rxn_rates(self.rxn_rates)
            
            
        #Step 2: simluate transition from parent steady state to daughter steady state
            #a. here we check when the system reaches a new steady state and sample before then
            #b. if steady states are similar mbe set simulation time?
            
        #try t-test for mean for all genes (for each thousand timepoints) compared to calculated steady state for current cell_type

        if len(init_mol_count) != 0:
            simulator.set_curr_mol_count(init_mol_count)

        f = open(count_csv, "w")
        f.truncate()
        f.close()      
        simulator.set_simulation_params(1000, True)
        simulator.simulate(list(np.random.randint(low = 1, high = 1000000, size = 4)))
        sim_out = pd.read_csv(count_csv)
        
        #identify molecules of interest in csv

        #steady state detection

        #collect Traisition samples
        if collect_samples:
            sim_out = pd.read_csv(count_csv)
            sim_out = sim_out.loc[:, ~sim_out.columns.str.contains('^Unnamed')]
            species = [col for col in sim_out.columns if col[-len(molecule):] == molecule]
            samples = np.random.random(size = num_samples) *1000
            time = list(np.cumsum(sim_out[sim_out.columns[0]]))
            time.append(np.inf)
            idxs = np.searchsorted(time, samples,side = 'right') - 1
            sample_out = sim_out.iloc[idxs][species]
            sample_out['Cell Type'] = [str(self.lineageName)+'T'] * num_samples
            sample_out.to_csv(sample_csv, mode = 'a', header = False)
        
        #Step 3: Sample num_sample transition cells and store to output and run sim_cell_type
        self.sim_cell_type(num_samples, simulator, molecule, count_csv, sample_csv)
        
        