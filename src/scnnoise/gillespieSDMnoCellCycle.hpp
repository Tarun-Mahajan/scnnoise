/* Copyright (C) 2021 Tarun Mahajan - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the GNU AGPLv3 license.
 *
 * You should have received a copy of the GNU AGPLv3 license with
 * this file. If not, please write to: tarunm3@illinois.edu, or 
 * visit : https://www.gnu.org/licenses/agpl-3.0.en.html
 */

// Optimized gillespie header file
#ifndef SDMnoCellCycle_H
#define SDMnoCellCycle_H

#include <vector>
#include <string>
#include "gillespieSSA.hpp"

namespace ScnnoiseInterface {
    /********************************************//**
    \brief A class for Gillespie's stochastic simulation algorithm.

    The Gillespie class creates an object to perform
    exact stochastic simulation for any given chemical
    reaction network.
    ***********************************************/

    class gillespieSDMnoCellCycle : public GillespieSSA {
    private:

    public:
    /* Member functions */
    // Constructor
        gillespieSDMnoCellCycle (int num_genes, std::string gene_filepath,
            std::string molecule_count_filepath,
            std::string count_save_file, bool keep_GRN,
            std::string GRN_filepath, int num_timepoints_save);

        void sort_reaction (int &rxn_selected) override;
        void update_cell_cycle_state (double next_time,
            double cur_time, RNG &generator) override;
        void init_cell_cycle_state (RNG &generator, double cur_time) override;

        std::string get_cur_cell_cycle_state () override;
    };
}

#endif
