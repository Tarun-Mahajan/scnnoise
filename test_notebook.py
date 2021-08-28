from scnnoise import _scnnoise
num_genes = 1
gene_filepath = "data/toy_data/single_gene_cell_cycle.csv"
molecule_count_filepath = "data/toy_data/single_gene_cell_cycle_molecule_count.csv"
count_save_file = "data/toy_data/single_gene_cell_cycle_count_output.csv"
keep_GRN = False
GRN_filepath = "dummy"
gene_obj = _scnnoise.gillespieSDMCellCycle(num_genes, gene_filepath, molecule_count_filepath, count_save_file, keep_GRN, GRN_filepath)
