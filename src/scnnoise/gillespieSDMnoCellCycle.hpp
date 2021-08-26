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
            std::string GRN_filepath);

        void sort_reaction (int &rxn_selected) override;
        void update_cell_cycle_state (double cur_time, RNG &generator) override;
        void init_cell_cycle_state (RNG &generator) override;
    };
}

#endif
