#include <iostream>
#include <string>
#include <vector>
#include "gillespieSDMCellCycle.hpp"
using namespace ScnnoiseInterface;

int main () {
    int num_genes = 1;
    std::string gene_filepath = "../data/toy_data/single_gene_cell_cycle.csv";
    std::string molecule_count_filepath = "../data/toy_data/single_gene_cell_cycle_molecule_count.csv";
    std::string count_save_file = "../data/toy_data/single_gene_cell_cycle_count_output.csv";
    bool keep_GRN = false;
    std::string  GRN_filepath = "dummy";
    double max_time = 20000;
    bool save_timeseries = true;
    int num_timepoints_save = 1000;
    gillespieSDMCellCycle obj(num_genes, gene_filepath,
        molecule_count_filepath, count_save_file, keep_GRN, GRN_filepath);

    std::vector<double> dosage_compensation(num_genes, 0.5);
    obj.set_cell_cycle_params();
    obj.set_dosage_compensation(dosage_compensation);
    obj.set_simulation_params (max_time, save_timeseries, num_timepoints_save);
    obj.simulate();
    // std::cout << obj.rxn_order[obj.rxn_order_map["A"]["gene off"]].propensity_val << std::endl;
    // std::cout << obj.reactions[0].propensity_vals["gene off"] << std::endl;
}
