#include <vector>
#include "gillespieSSA.hpp"
#include "gillespieSDMnoCellCycle.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <random>
#include <math.h>
// #include "graph.hpp"
#include "scnnoise.hpp"


namespace ScnnoiseInterface {
    /* function definitions */
    // Constructor
    gillespieSDMnoCellCycle::gillespieSDMnoCellCycle (int num_genes,
        std::string gene_filepath,
        std::string molecule_count_filepath,
        std::string count_save_file, bool keep_GRN,
        std::string GRN_filepath):
    GillespieSSA (num_genes, gene_filepath,
        molecule_count_filepath,
        count_save_file, keep_GRN,
        GRN_filepath) {
    }

    void gillespieSDMnoCellCycle::sort_reaction (int &rxn_selected) {
        if (rxn_selected > 0) {
            std::string gene_name = gene_map[rxn_order[rxn_selected].gene_id];
            std::string rxn_name = rxn_order[rxn_selected].rxn_name;
            rxn_order_map[gene_name][rxn_name] -= 1;

            gene_name = gene_map[rxn_order[rxn_selected - 1].gene_id];
            rxn_name = rxn_order[rxn_selected - 1].rxn_name;
            rxn_order_map[gene_name][rxn_name] += 1;

            std::swap(rxn_order[rxn_selected - 1], rxn_order[rxn_selected]);
            rxn_selected -= 1;
        }
    }

    inline void gillespieSDMnoCellCycle::update_cell_cycle_state (double next_time,
        double cur_time, RNG &generator) {

    }

    void gillespieSDMnoCellCycle::init_cell_cycle_state (RNG &generator, double cur_time) {

    }
}
