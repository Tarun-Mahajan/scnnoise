#include <iostream>
#include <string>
#include <vector>
#include <random>
#include "gillespieSDMCellCycle.hpp"
#include "gillespieSDMnoCellCycle.hpp"
using namespace ScnnoiseInterface;

int main () {
    std::random_device rd;
    // std::vector<std::uint_least32_t> rd_seeds = {rd(), rd(), rd(), rd()};
    std::vector<unsigned int> random_seeds {582654328, 1065236345, 322147403, 2229968939};
    std::vector<std::uint_least32_t> rd_seeds =
        {582654328, 1065236345, 322147403, 2229968939};
    std::seed_seq sd(rd_seeds.begin(), rd_seeds.end());
    thread_local static std::mt19937 generator{sd};
    std::discrete_distribution<int> distribution {0.3, 0.7};

    // std::cout << rd() << std::endl;

    int num_genes = 2;
    /*
    std::string gene_filepath = "../data/toy_data/single_gene_cell_cycle.csv";
    std::string molecule_count_filepath = "../data/toy_data/single_gene_cell_cycle_molecule_count.csv";
    std::string count_save_file_before =
        "../data/output/cell_cycle_arrest_test/single_gene_cell_cycle_count_output_before_arrest";
    std::string count_save_file_after =
        "../data/output/cell_cycle_arrest_test/single_gene_cell_cycle_count_output_after_arrest";
    std::string count_save_file;
    std::string connector_ = "_";
    std::string file_extension = ".csv";
    bool keep_GRN = false;
    std::string  GRN_filepath = "dummy";
    double max_time = 20000;
    bool save_timeseries = true;
    int num_timepoints_save = 1000;
    unsigned int nsamples = 1000;

    for (unsigned int b = 0; b < nsamples; ++b) {
        std::cout << "sample no. = " << b << std::endl;
        count_save_file = count_save_file_before +
            connector_ + std::to_string(b) + file_extension;

        gillespieSDMCellCycle obj(num_genes, gene_filepath,
            molecule_count_filepath, count_save_file, keep_GRN, GRN_filepath);

        std::vector<double> dosage_compensation(num_genes, 0.7);
        obj.set_cell_cycle_params();
        obj.set_dosage_compensation(dosage_compensation);
        obj.set_simulation_params (max_time, save_timeseries, num_timepoints_save);
        obj.set_cell_cycle_frozen_state();
        obj.simulate();

        count_save_file = count_save_file_after +
            connector_ + std::to_string(b) + file_extension;
        obj.change_output_filepath(count_save_file);
        int arrest_  = distribution(generator);
        if (arrest_ == 1) {
            obj.set_cell_cycle_frozen_state(true, "G2");
        }else{
            obj.set_cell_cycle_frozen_state();
        }

        obj.simulate();
    }**/

    /*
    // single gene two-state reduced model cell-cycle
    std::string gene_filepath = "../data/input/single_gene_reduced_model_cell_cycle/single_gene_reduced_model_cell_cycle.csv";
    std::string molecule_count_filepath = "../data/input/single_gene_reduced_model_cell_cycle/single_gene_reduced_model_cell_cycle_molecule_count.csv";
    std::string burst_size_filepath = "../data/input/single_gene_reduced_model_cell_cycle/single_gene_reduced_model_cell_cycle_burst_size.csv";
    std::string count_save_file_before =
        "../data/output/single_cell_reduced_model_cell_cycle/single_gene_reduced_model_cell_cycle_count_output_before_arrest";
    std::string count_save_file_after =
        "../data/output/single_cell_reduced_model_cell_cycle/single_gene_reduced_model_cell_cycle_count_output_after_arrest";
    std::string count_save_file;
    std::string connector_ = "_";
    std::string file_extension = ".csv";
    bool keep_GRN = false;
    std::string  GRN_filepath = "dummy";
    double max_time = 20000;
    bool save_timeseries = true;
    int num_timepoints_save = 1000;
    unsigned int nsamples = 1000;

    for (unsigned int b = 0; b < nsamples; ++b) {
        std::cout << "sample no. = " << b << std::endl;
        count_save_file = count_save_file_before +
            connector_ + std::to_string(b) + file_extension;

        gillespieSDMCellCycle obj(num_genes, gene_filepath,
            molecule_count_filepath, count_save_file, keep_GRN, GRN_filepath);

        std::vector<double> dosage_compensation(num_genes, 0.5);
        obj.set_reduced_model_stoichio_factor(burst_size_filepath);
        obj.set_cell_cycle_params();
        obj.set_dosage_compensation(dosage_compensation);
        obj.set_simulation_params (max_time, save_timeseries, num_timepoints_save);
        obj.set_cell_cycle_frozen_state();
        obj.simulate();

        count_save_file = count_save_file_after +
            connector_ + std::to_string(b) + file_extension;
        gillespieSDMCellCycle obj1(num_genes, gene_filepath,
            molecule_count_filepath, count_save_file, keep_GRN, GRN_filepath);
        // obj.change_output_filepath(count_save_file);
        // int arrest_  = distribution(generator);
        // if (arrest_ == 1) {
        //     obj.set_cell_cycle_frozen_state(true, "G2");
        // }else{
        //     obj.set_cell_cycle_frozen_state();
        // }
        obj1.set_reduced_model_stoichio_factor(burst_size_filepath);
        obj1.set_cell_cycle_params();
        obj1.set_dosage_compensation(dosage_compensation);
        obj1.set_simulation_params (max_time, save_timeseries, num_timepoints_save);
        obj1.set_cell_cycle_frozen_state(true, "G2");
        obj1.simulate();
    }**/

    // two-gene cascade two-state reduced model without cell-cycle
    std::string gene_filepath = "../data/input/two_gene_cascade/two_gene_cascade_gene_state_two_state_reduced.csv";
    std::string molecule_count_filepath = "../data/input/two_gene_cascade/two_gene_cascade_two_state_reduced_molecule_count.csv";
    std::string burst_size_filepath = "../data/input/two_gene_cascade/two_gene_cascade_two_state_reduced_burst_size.csv";
    std::string count_save_file_before =
        "../data/output/two_gene_cascade/two_gene_cascade_two_state_reduced_count_output";
    std::string count_save_file_after =
        "../data/output/two_gene_cascade/single_gene_two_state_cell_cycle_count_output_after_treatment";
    std::string count_save_file;
    std::string connector_ = "_";
    std::string file_extension = ".csv";
    bool keep_GRN = true;
    std::string  GRN_filepath = "../data/input/two_gene_cascade/two_gene_cascade_GRN.csv";
    double max_time = 10000;
    bool save_timeseries = true;
    int num_timepoints_save = 10000;
    unsigned int nsamples = 1;

    for (unsigned int b = 0; b < nsamples; ++b) {
        std::cout << "sample no. = " << b << std::endl;
        count_save_file = count_save_file_before +
            connector_ + std::to_string(b) + file_extension;

        gillespieSDMCellCycle obj(num_genes, gene_filepath,
            molecule_count_filepath, count_save_file, keep_GRN, GRN_filepath, num_timepoints_save);
        std::vector<double> dosage_compensation(num_genes, 0.5);
        obj.set_reduced_model_stoichio_factor(burst_size_filepath);
        obj.set_cell_cycle_params();
        obj.set_dosage_compensation(dosage_compensation);
        obj.set_simulation_params (max_time, save_timeseries);
        obj.set_cell_cycle_frozen_state("true", "G1");
        obj.simulate(random_seeds);

        // After treatment
        // count_save_file = count_save_file_after +
        //     connector_ + std::to_string(b) + file_extension;
        // obj.change_output_filepath(count_save_file);
        // // obj.set_cell_cycle_params("timing", "gamma", distribution_params_init,
        // //     0.2);
        // obj.set_cell_cycle_frozen_state(true, "G2");
        // obj.simulate();


    }
    // gillespieSDMCellCycle obj(num_genes, gene_filepath,
    //     molecule_count_filepath, count_save_file, keep_GRN, GRN_filepath);
    //
    // std::vector<double> dosage_compensation(num_genes, 0.5);
    // obj.set_cell_cycle_params();
    // obj.set_dosage_compensation(dosage_compensation);
    // obj.set_simulation_params (max_time, save_timeseries, num_timepoints_save);
    // obj.set_cell_cycle_frozen_state();
    // obj.simulate();
    //
    // count_save_file = "../data/output/cell_cycle_arrest_test/single_gene_cell_cycle_count_output_after_arrest.csv";
    // obj.change_output_filepath(count_save_file);
    // obj.set_cell_cycle_frozen_state(true, "G2");
    // obj.simulate();
    // std::cout << obj.rxn_order[obj.rxn_order_map["A"]["gene off"]].propensity_val << std::endl;
    // std::cout << obj.reactions[0].propensity_vals["gene off"] << std::endl;
}
