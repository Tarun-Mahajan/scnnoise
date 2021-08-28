#include <vector>
#include "gillespieSSA.hpp"
#include "gillespieSDMCellCycle.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <random>
#include <math.h>
// #include "graph.hpp"
// #include "scnnoise.hpp"


namespace ScnnoiseInterface {
    /* function definitions */
    // Constructor
    gillespieSDMCellCycle::gillespieSDMCellCycle (int num_genes, std::string gene_filepath,
        std::string molecule_count_filepath,
        std::string count_save_file, bool keep_GRN,
        std::string GRN_filepath):
    GillespieSSA (num_genes, gene_filepath,
        molecule_count_filepath,
        count_save_file, keep_GRN,
        GRN_filepath) {
    }

    void gillespieSDMCellCycle::swap_rxn_order (rxn_order_struct &A,
        rxn_order_struct &B) {
        rxn_order_struct temp = A;
        A = B;
        B = temp;
    }

    void gillespieSDMCellCycle::sort_reaction (int &rxn_selected) {
        if (rxn_selected > 0) {
            std::string gene_name = gene_map[rxn_order[rxn_selected].gene_id];
            std::string rxn_name = rxn_order[rxn_selected].rxn_name;
            rxn_order_map[gene_name][rxn_name] -= 1;

            gene_name = gene_map[rxn_order[rxn_selected - 1].gene_id];
            rxn_name = rxn_order[rxn_selected - 1].rxn_name;
            rxn_order_map[gene_name][rxn_name] += 1;

            swap_rxn_order(rxn_order[rxn_selected - 1], rxn_order[rxn_selected]);
            rxn_selected -= 1;
        }
    }

    void gillespieSDMCellCycle::update_propensity_cell_cycle () {
        for (std::size_t gene = 0; gene < num_genes; ++gene) {
            std::string gene_type = reactions[gene].gene_type;
            auto it_on = gene_type_info[gene_type].species_rev_map.find("gene on");
            if (it_on != gene_type_info[gene_type].species_rev_map.end()) {
                // Update propensity for gene on reactions
                total_propensity -= reactions[gene].propensity_vals["gene on"];
                double new_propensity =
                    compute_propensity(gene_map[gene], "gene on");
                reactions[gene].propensity_vals["gene on"] = new_propensity;
                unsigned int rxn_order_id =
                    rxn_order_map[gene_map[gene]]["gene on"];
                rxn_order[rxn_order_id].propensity_val = new_propensity;
                total_propensity += new_propensity;

                // Update propensity for gene off reactions
                total_propensity -= reactions[gene].propensity_vals["gene off"];
                new_propensity =
                    compute_propensity(gene_map[gene], "gene off");
                reactions[gene].propensity_vals["gene off"] = new_propensity;
                rxn_order_id =
                    rxn_order_map[gene_map[gene]]["gene off"];
                rxn_order[rxn_order_id].propensity_val = new_propensity;
                total_propensity += new_propensity;
            }else {
            }
            // Update propensity for transcription reactions
            total_propensity -= reactions[gene].propensity_vals["transcription"];
            double new_propensity =
                compute_propensity(gene_map[gene], "transcription");
            reactions[gene].propensity_vals["transcription"] = new_propensity;
            unsigned int rxn_order_id =
                rxn_order_map[gene_map[gene]]["transcription"];
            rxn_order[rxn_order_id].propensity_val = new_propensity;
            total_propensity += new_propensity;
        }
    }

    void gillespieSDMCellCycle::replicate_genes () {
        for (std::size_t gene = 0; gene < num_genes; ++gene) {
            std::string gene_type = reactions[gene].gene_type;
            auto it_on = gene_type_info[gene_type].species_rev_map.find("gene on");
            if (it_on != gene_type_info[gene_type].species_rev_map.end()) {
                int on_id = gene_type_info[gene_type].species_rev_map["gene on"];
                int off_id = gene_type_info[gene_type].species_rev_map["gene off"];
                int gene_copy_number =
                    reactions[gene].molecule_count_cur[on_id] +
                    reactions[gene].molecule_count_cur[off_id];
                reactions[gene].molecule_count_cur[off_id] = 2*gene_copy_number;
                reactions[gene].molecule_count_cur[on_id] = 0;
            }else {
                std::string search_string = "two-state reduced";
                size_t found = gene_type.find(search_string);
                if (found != std::string::npos) {
                    double &stoichio_ =
                        stoichio_factors[gene].rxns["transcription"].products_factors["mRNA"];
                    stoichio_ = 2;
                }else{
                    search_string = "constitutive";
                    found = gene_type.find(search_string);
                    if (found != std::string::npos) {
                        reactions[gene].rxn_rates["transcription"] *= 2;
                    }
                }
            }
        }
    }

    void gillespieSDMCellCycle::cell_division (RNG& generator) {
        for (std::size_t gene = 0; gene < num_genes; ++gene) {
            std::string gene_type = reactions[gene].gene_type;
            auto it_on = gene_type_info[gene_type].species_rev_map.find("gene on");
            int on_id = -1;
            int off_id = -1;
            if (it_on != gene_type_info[gene_type].species_rev_map.end()) {
                on_id = gene_type_info[gene_type].species_rev_map["gene on"];
                off_id = gene_type_info[gene_type].species_rev_map["gene off"];
                int gene_copy_number =
                    reactions[gene].molecule_count_cur[on_id] +
                    reactions[gene].molecule_count_cur[off_id];
                gene_copy_number = gene_copy_number/2;
                reactions[gene].molecule_count_cur[off_id] =
                    floor(reactions[gene].molecule_count_cur[off_id]/2);
                reactions[gene].molecule_count_cur[on_id] =
                    floor(reactions[gene].molecule_count_cur[on_id]/2);
                int copy_number_diff = gene_copy_number - (
                    reactions[gene].molecule_count_cur[off_id] +
                    reactions[gene].molecule_count_cur[on_id]);
                if (copy_number_diff > 0) {
                    std::uniform_int_distribution<int> distribution(0,1);
                    for (int i = 0; i < copy_number_diff; ++i) {
                        if (distribution(generator) == 0) {
                            reactions[gene].molecule_count_cur[off_id] += 1;
                        }else{
                            reactions[gene].molecule_count_cur[on_id] += 1;
                        }
                    }
                }
            }else {
                std::string search_string = "two-state reduced";
                std::size_t found = gene_type.find(search_string);
                if (found != std::string::npos) {
                    double &stoichio_ =
                        stoichio_factors[gene].rxns["transcription"].products_factors["mRNA"];
                    stoichio_ = 1;
                }else{
                    search_string = "constitutive";
                    found = gene_type.find(search_string);
                    if (found != std::string::npos) {
                        reactions[gene].rxn_rates["transcription"] /= 2;
                    }
                }
            }

            for (std::size_t species = 0;
                species < reactions[gene].molecule_count_cur.size(); ++species) {
                if (species != off_id && species != on_id) {
                    int molecule_count = reactions[gene].molecule_count_cur[species];
                    reactions[gene].molecule_count_cur[species] = 0;
                    std::uniform_int_distribution<int> distribution(0,1);
                    for (int mol_count = 0; mol_count < molecule_count; ++mol_count) {
                        if (distribution(generator) == 1) {
                            reactions[gene].molecule_count_cur[species] += 1;
                        }
                    }
                }
            }
        }
    }

    void gillespieSDMCellCycle::perform_dosage_compensation () {
        for (std::size_t gene = 0; gene < num_genes; ++gene) {
            std::string gene_type = reactions[gene].gene_type;
            auto it_on = gene_type_info[gene_type].species_rev_map.find("gene on");
            if (it_on != gene_type_info[gene_type].species_rev_map.end()) {
                reactions[gene].rxn_rates["gene on"] *= dosage_compensation[gene];
            }else {
                std::string search_string = "two-state reduced";
                size_t found = gene_type.find(search_string);
                if (found != std::string::npos) {
                    reactions[gene].rxn_rates["transcription"] *=
                        dosage_compensation[gene];
                }
            }
        }
    }

    void gillespieSDMCellCycle::remove_dosage_compensation () {
        for (std::size_t gene = 0; gene < num_genes; ++gene) {
            std::string gene_type = reactions[gene].gene_type;
            auto it_on = gene_type_info[gene_type].species_rev_map.find("gene on");
            if (it_on != gene_type_info[gene_type].species_rev_map.end()) {
                reactions[gene].rxn_rates["gene on"] /= dosage_compensation[gene];
            }else {
                std::string search_string = "two-state reduced";
                size_t found = gene_type.find(search_string);
                if (found != std::string::npos) {
                    reactions[gene].rxn_rates["transcription"] /=
                        dosage_compensation[gene];
                }
            }
        }
    }

    void gillespieSDMCellCycle::check_if_replication () {
        if (cur_time < cell_cycle_start_time +
            cell_cycle_length * cell_cycle_params.replication_time_factor &&
            next_time > cell_cycle_start_time +
            cell_cycle_length * cell_cycle_params.replication_time_factor) {
                replicate_genes();
                perform_dosage_compensation();
                update_propensity_cell_cycle();
                // compute_total_propensity();
                current_cell_cycle_state = "G2";
            }
    }

    void gillespieSDMCellCycle::check_if_division (RNG& generator) {
        if (cur_time < cell_cycle_start_time + cell_cycle_length &&
            next_time > cell_cycle_start_time + cell_cycle_length) {
                cell_division(generator);
                remove_dosage_compensation();
                update_propensity_cell_cycle();
                // compute_total_propensity();
                cell_cycle_start_time = cell_cycle_start_time + cell_cycle_length;
                sample_cell_cycle_time(generator);
                current_cell_cycle_state = "G1";
            }
    }

    void gillespieSDMCellCycle::update_cell_cycle_state (double next_time,
        double cur_time,
        RNG& generator) {
        this->next_time = next_time;
        this->cur_time = cur_time;
        if (!is_frozen_cell_cycle) {
            check_if_replication();
            check_if_division(generator);
        }
    }

    void gillespieSDMCellCycle::set_cell_cycle_params (std::string mechanism,
        std::string distribution_,
        std::map<std::string, double> distribution_params,
        double replication_time_factor) {
        cell_cycle_params.mechanism = mechanism;
        cell_cycle_params.distribution_ = distribution_;
        cell_cycle_params.distribution_params = distribution_params;
        cell_cycle_params.replication_time_factor = replication_time_factor;
    }

    void gillespieSDMCellCycle::set_dosage_compensation (
        std::vector<double> dosage_compensation) {
        this->dosage_compensation = dosage_compensation;
    }

    void gillespieSDMCellCycle::sample_cell_cycle_time (RNG &generator) {
        cell_cycle_length = 0;
        if (cell_cycle_params.mechanism == "timing") {
            if (cell_cycle_params.distribution_ == "gamma") {
                double alpha_ = cell_cycle_params.distribution_params["alpha"];
                double beta_ = cell_cycle_params.distribution_params["beta"];
                thread_local std::gamma_distribution<double> distribution(alpha_,
                    beta_);
                while (cell_cycle_length <= 0) {
                    cell_cycle_length = distribution(generator);
                }
            }
        }
    }

    void gillespieSDMCellCycle::set_cur_time (double cur_time) {
        this->cur_time = cur_time;
    }

    void gillespieSDMCellCycle::set_cell_cycle_start_time () {
        this->cell_cycle_start_time = this->cur_time;
    }

    void gillespieSDMCellCycle::move_to_G1 (RNG &generator) {
        cell_division(generator);
        remove_dosage_compensation();
        update_propensity_cell_cycle();
        compute_total_propensity();
        current_cell_cycle_state = "G1";
    }

    void gillespieSDMCellCycle::move_to_G2 () {
        replicate_genes();
        perform_dosage_compensation();
        update_propensity_cell_cycle();
        compute_total_propensity();
        current_cell_cycle_state = "G2";
    }

    void gillespieSDMCellCycle::init_cell_cycle_state (RNG &generator,
        double cur_time) {
        if (!is_frozen_cell_cycle) {
            if (current_cell_cycle_state.empty()) {
                current_cell_cycle_state = "G1";
            }else{
                if (current_cell_cycle_state == "G2") {
                    move_to_G1(generator);
                }
            }
        }else{
            if (current_cell_cycle_state.empty()) {
                if (frozen_state == "G1") {
                    current_cell_cycle_state = "G1";
                }else{
                    move_to_G2();
                }
            }else{
                if (current_cell_cycle_state == "G1") {
                    if (frozen_state == "G2") {
                        move_to_G2();
                    }
                }else{
                    if (frozen_state == "G1") {
                        move_to_G1(generator);
                    }
                }
            }
        }
        sample_cell_cycle_time (generator);
        set_cur_time (cur_time);
        set_cell_cycle_start_time ();
    }

    void gillespieSDMCellCycle::set_cell_cycle_frozen_state (bool is_frozen_cell_cycle,
        std::string frozen_state) {
        this->is_frozen_cell_cycle = is_frozen_cell_cycle;
        this->frozen_state = frozen_state;
    }
}
