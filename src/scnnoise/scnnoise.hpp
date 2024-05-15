// scNNoiSE header file
#ifndef SCNNOISE
#define SCNNOISE

#include "graph.hpp"
#include "graph_derived.hpp"
#include <vector>
#include <map>
#include <random>
#include <string>
#include <math.h>
// #include "graph.hpp"



/**
 * \file scnnoise.hpp
 * \brief Contains the definition of various structs and classes used in the scNNoiSE library.
 */


namespace ScnnoiseInterface {

    /**
     * \struct rxn_order_struct
     * \brief A struct for reaction order.
     *
     * This struct represents the order of a reaction and contains the following members:
     * - gene_id: Numeric ID of the gene.
     * - rxn_name: Name of the reaction.
     * - propensity_val: Propensity value of the reaction.
     */
    struct rxn_order_struct {
        int gene_id;
        std::string rxn_name;
        double propensity_val;
    };

    /**
     * \struct rxn_struct
     * \brief A struct for a chemical reaction.
     *
     * This struct represents a chemical reaction and contains the following members:
     * - reactants_stoichio: A map of reactant species and their stoichiometric coefficients.
     * - products_stoichio: A map of product species and their stoichiometric coefficients.
     */
    struct rxn_struct {
        // std::vector<int> reactants;
        // std::vector<int> products;
        std::map<std::string, int> reactants_stoichio;
        std::map<std::string, int> products_stoichio;
        // double rxn_rate;
        // double propensity_val;
    };

    /**
     * \struct gene_type_struct
     * \brief A struct representing a gene type.
     *
     * This struct represents a gene type and contains the following members:
     * - num_rxns: Number of reactions associated with the gene type.
     * - num_species: Number of species associated with the gene type.
     * - species_rev_map: A map of species names to their integer IDs.
     * - species_map: A map of integer IDs to species names.
     * - rxn_rev_map: A map of reaction names to their integer IDs.
     * - rxn_map: A map of integer IDs to reaction names.
     * - rxns: A map of reaction names to their corresponding rxn_struct.
     * - gene_rxn_dependency: A vector of GraphSpace::GraphDependency representing the dependency graph of gene reactions.
     */
    struct gene_type_struct {
        unsigned int num_rxns;
        unsigned int num_species;
        std::map<std::string, int> species_rev_map;
        std::map<int, std::string> species_map;
        std::map<std::string, int> rxn_rev_map;
        std::map<int, std::string> rxn_map;
        std::map<std::string, rxn_struct> rxns;
        std::vector<GraphSpace::GraphDependency> gene_rxn_dependency;
    };

    /**
     * \struct stoichio_factor_species_struct
     * \brief A struct representing stoichiometric factors for a species in a reaction.
     *
     * This struct represents the stoichiometric factors for a species in a reaction and contains the following members:
     * - reactants_factors: A map of reactant species and their stoichiometric factors.
     * - products_factors: A map of product species and their stoichiometric factors.
     */
    struct stoichio_factor_species_struct {
        // std::vector<int> reactants;
        // std::vector<int> products;
        std::map<std::string, double> reactants_factors;
        std::map<std::string, double> products_factors;
        // double rxn_rate;
        // double propensity_val;
    };

    /**
     * \struct stoichio_factor_struct
     * \brief A struct representing stoichiometric factors for reactions.
     *
     * This struct represents the stoichiometric factors for reactions and contains the following members:
     * - rxns: A map of reaction names to their corresponding stoichio_factor_species_struct.
     */
    struct stoichio_factor_struct {
        // std::vector<int> reactants;
        // std::vector<int> products;
        std::map<std::string, stoichio_factor_species_struct> rxns;
        // double rxn_rate;
        // double propensity_val;
    };

    /**
     * \struct gene_rxn_channel_struct
     * \brief A struct representing information for all reaction channels for a gene.
     *
     * This struct represents information for all reaction channels for a gene and contains the following members:
     * - gene_name: Name of the gene.
     * - gene_type: Type of the gene.
     * - GRN_rxn_IN: A vector of reaction names for which reaction rate is affected by the GRN.
     * - GRN_species_OUT: A vector of species names which cause regulation in the GRN.
     * - rxn_rates: A map of reaction names to their corresponding rates.
     * - molecule_count_cur: A vector of current molecule counts for chemical species related to the gene.
     * - propensity_vals: A map of reaction names to their corresponding propensity values.
     */
    struct gene_rxn_channel_struct {
        std::string gene_name;
        std::string gene_type;
        std::vector<std::string> GRN_rxn_IN;
        std::vector<std::string> GRN_species_OUT;
        // std::vector<std::string> rxn_names;
        // std::vector<int> rxn_rates;
        std::map<std::string, double> rxn_rates;
        std::vector<int> molecule_count_cur;
        std::map<std::string, double> propensity_vals;
    };

    /**
     * \struct molecule_history_struct
     * \brief A struct representing the history of molecule counts.
     *
     * This struct represents the history of molecule counts and contains the following members:
     * - molecule_count: A vector of vectors representing the molecule counts at different time points.
     */
    struct molecule_history_struct {
        std::vector<std::vector<int>> molecule_count;
    };

    class scNNoiSE {
    public:
        typedef std::mt19937 RNG; // Mersenne Twister 19937 generator
        std::vector<RNG> generator; // random number generator
        int num_rxns; // number of reactions
        int num_genes; // number of genes
        std::map<std::string, gene_type_struct> gene_type_info; // gene type information
        std::map<int, std::string> gene_map; // gene map
        std::map<std::string, int> gene_rev_map; // gene reverse map

        /********************************************//**
         \brief Dependency graph for reaction channels for
                different combinations of gene types and num
                of alternatively spliced mRNA.
         ***********************************************/
        // std::vector<GraphSpace::GraphDependency> gene_rxn_dependency;

        /********************************************//**
         \brief Number of chemical species for each
                different combinations of gene types and num
                of alternatively spliced mRNA.
         ***********************************************/
        std::vector<int> num_species_gene_type;

        /********************************************//**
        \brief Number of chemical reactions for each
             different combinations of gene types and num
             of alternatively spliced mRNA.
        ***********************************************/
        std::vector<int> num_rxns_gene_type;

        /********************************************//**
        \brief Reaction search order

        Order in which reaction propensities are added
        for reaction channel selection at each step of the
        simulation.
        ***********************************************/
        std::vector<rxn_order_struct> rxn_order;

        std::map<std::string, std::map<std::string, unsigned int>> rxn_order_map;

        /********************************************//**
        \brief Gene regulatory network.

        Gene regulatory network is stored in a vector of
        type graph, which is a class for gene regulatory
        networks.
        ***********************************************/
        std::vector<GraphSpace::GRN> network;

        bool keep_GRN;
        
        std::vector<gene_rxn_channel_struct> reactions;

        std::vector<stoichio_factor_struct> stoichio_factors;

        std::vector<std::string> burst_size_distribution;

        std::vector<double> gene_burst_sizes;

        std::vector<double> burst_sizes_mean;

        std::vector<double> gene_copy_number;

        std::map<int, std::map<std::string, double>> max_rxn_rate_change;

        std::vector<double> basal_regulation;

        std::vector<double> max_freq_regulation;

        std::string regulation_type;

        // data members for steady state testing of moments
        bool is_steady_state_reached;

        unsigned int num_history_statistics;

        unsigned int burn_in_rxn_count;

        unsigned int after_burn_in_rxn_count;

        std::vector<std::vector<std::vector<double>>> history_statistics;


        // /********************************************//**
        //  \brief Vector for status in GRN
        //
        //  Vector to store status of the molecular species--
        //  whether the species participates in the GRN or not.
        //  ***********************************************/
        // std::vector<bool> in_GRN;

        /********************************************//**
        \brief Simulator settings.

        Gene regulatory network is stored in a vector of
        type GRN, which is a class for gene regulatory
        networks.
        ***********************************************/
        // std::vector<Simulator_Preprocess> simulator_settings;

        /********************************************//**
        \brief Simulator settings.

        Gene regulatory network is stored in a vector of
        type GRN, which is a class for gene regulatory
        networks.
        ***********************************************/
        // std::vector<Simulator> simulation;

        // Gillespie simulator;

        /********************************************//**
        \brief Time period for simulation run.

        The total time period for which the chemical
        reaction network should be simulated.
        ***********************************************/
        double max_time;

        double burn_in;

        double num_points_to_collect;

        double time_interval_to_save;

        bool save_at_time_interval;

        bool save_at_random_times;

        std::vector<double> random_times_to_save;

        /********************************************//**
        \brief Save time series count data.

        This boolean variable decides whether time series
        count data should be stored or not. Default value
        is false, when only steady state count values are
        reported.
        ***********************************************/
        bool save_timeseries;

        bool save_timeseries_all;

        /********************************************//**
        \brief Number of previous time points for which molecular counts are stored in
            a vector of vectors..

        If save_timeseries is set to True, num_timepoints_save gives the number of
        past values of molecular counts to store in an 2D array
        (vector of vectors) of size . Once num_timepoints_save past values have been
        stored for each moleulcar species, the contents of the array are
        appened to an output file, and the array is emptied. After this,
        array starts storing the futures values of molecular counts until
        num_timepoints_save are stored, when the values are again appened to the
        output file. This process continues until the simulation terminates.
        At termination, the array again appends its values to the output file.
        ***********************************************/
        int num_timepoints_save;

        std::vector<std::map<std::string, unsigned int>> count_rxns_fired;

        bool count_rxns;
        unsigned int stop_rxn_count;

        double total_propensity;

        /********************************************//**
        \brief Vector of vectors to store running molecule counts for all species
            for last 'num_timepoints_save' time points.

        A 2D array stored as a vector of vectors to store running molecular counts
        for all the molecular specis in the system. This array ensures that the
        molecular counts are written to an output file at a fixed interval. If
        num_timepoints_save > 1, This reduces the number of times data has to be
        written to the output file.
        ***********************************************/
        std::vector<std::vector<std::vector<int>>> molecule_count_history;

        std::vector<std::string> cell_cycle_phase_history;

        /********************************************//**
        \brief Vector to store time points along the simulation path.
        ***********************************************/
        std::vector<double> time_history;

        std::vector<std::string> rxn_history;

        // file to save molecule count history
        std::string count_save_file;

        std::vector<std::vector<double>> running_mean;

        std::vector<std::vector<double>> running_var;

        std::vector<double> running_cov;

        // running marginal and joint probabilities
        std::vector<std::vector<std::map<unsigned int, double>>> running_marginal_probs;

        std::vector<std::map<std::string, double>> running_joint_probs;

        unsigned int running_probs_buffer_size;

        std::string filepath_probs;

        double total_time_norm;

        // public:
        /********************************************//**
        \brief Constructor for scNNoiSE.

        \param[in] num_rxns Number of channels in the
                list of chemical reaction channels.
        \param[in] num_species number of molecular species
                in the system.
        \param[in] num_species_gene_type vector containing
                the number of chemical species for each
                gene type.
        \param[in] num_rxns_gene_type vector containing
                the number of chemical reactions for each
                gene type.
        ***********************************************/
        scNNoiSE (int num_genes, std::string gene_filepath,
            std::string molecule_count_filepath,
            std::string count_save_file, bool keep_GRN,
            std::string GRN_filepath, int num_timepoints_save);

    /********************************************//**
     \brief Add state for a gene.

     Function to add state for a gene. Gene state includes
     gene-type, gene copy number and the number of unique alternatively
     spliced (AS) mRNA.
     \param[in] gene_id integer id for the gene.
     \param[in] gene_type type for the gene.
     \param[in] gene_copy_number Gene copy number.
     \param[in] num_splice_variants number of AS variants for the gene identified by
                gene_id.
     ***********************************************/

        std::string match_and_return_gene_type (std::string in_gene_type);

        void init_rxn_order ();

        void init_molecule_count (std::string filepath);

        void init_molecule_count_history ();

        void init_gene_states_from_file (std::string filepath);

        void create_init_gene_type_info ();

        gene_type_struct create_constitutive_type ();

        gene_type_struct create_constitutive_nascent_type ();

        gene_type_struct create_two_state_type ();

        gene_type_struct create_two_state_two_gene_cascade_activation_type ();

        gene_type_struct create_two_state_mRNA_type ();

        gene_type_struct create_two_state_nascent_type ();

        void set_reduced_model_stoichio_factor (std::string filepath);

        void set_reduced_model_burst_size_manual (int gene_id, double burst_size,
            double copy_number, std::string distribution_name);

        gene_type_struct create_two_state_reduced_type ();

        gene_type_struct create_two_state_reduced_nascent_type();

        gene_type_struct create_two_state_reduced_mRNA_type ();

        gene_type_struct create_two_state_reduced_nascent_mRNA_type ();

        void init_max_rxn_rate_change ();

        void create_GRN (std::string filepath);

        void create_GRN_from_file (std::string filepath);

        /********************************************//**
        \brief Add dependency edge.

        For any gene in the chemical reaction network, the
        dependency graph gives the directed relationships
        between the reactions for that gene corresponding
        to the type for that gene.

        \param[in] gene_type integer id for gene type
        \param[in] src source rxn for the edge.
        \param[in] dest destination rxn for the edge.
        ***********************************************/
        // void add_dependency_edge (int gene_type, int src, int dest);
        typedef std::map<std::string, std::map<std::string, int>> reactant_product_type;
        void add_new_dependency_graph (std::string gene_type_name,
          std::map<std::string, int> species_map, std::map<int, std::string> rxn_map,
          reactant_product_type rxns_reactants,
          reactant_product_type rxns_products, std::vector<std::vector<int>> edge_list);

        int factorial (int num);

        int factorial_ratio_propensity_func (int N, int r);

        void compute_total_propensity ();

        double compute_regulation_function (int gene_id, std::string rxn_name);

        double compute_propensity (std::string gene_name, std::string rxn_name);

        /********************************************//**
        \brief Simulating stochastic gene expression.

        A pure virtual function for simulating stochastic gene
        expression. Needs to be overridden in any derived class.
        ***********************************************/
        virtual void simulate (bool compute_statistics = false,
            std::string statistics_file = "dummy", bool verbose = true, 
            bool cell_cycle_sim_frozen=true, 
            bool init_dosage_comp_adj=true) = 0;

        void set_simulation_params (double max_time = 10000,
            bool save_timeseries = false);

        /********************************************//**
        \brief Function to compute gene expression regulation by transcription factors
        ***********************************************/
        // double regulation_function (int gene_selected, int rxn);

        /********************************************//**
        \brief Hill function for regulation
        ***********************************************/
        double hill_function (int tf_count, double hill_coeff, double half_maximal,
            bool activator, double prob_contr);

        void init_stoichio_factors();

        void change_output_filepath (std::string new_filepath);

        void swap_rxn_rates (std::map<std::string, std::map<std::string, double>> rxn_rates);

        void set_curr_mol_count (std::map<std::string, std::map<std::string, int>> init_count);

        void update_burst_size (RNG &generator, int rxn_selected);

        void set_save_timeseries_all (bool save_timeseries_all);

        void set_count_rxns_fired (bool count_rxns, bool compute_statistics,
            unsigned int stop_rxn_count = pow(10, 6));

        bool update_rxn_count (int rxn_selected, bool &stop_sim,
            bool &reached_rxn_count, bool compute_statistics,
            unsigned int burn_in_rxn_count);

        void update_burst_size_init ();

        void set_regulation_type (std::string regulation_type = "hill additive");

        void find_random_times_to_save (RNG &generator, double burn_in,
            double max_time);

        void set_num_points_to_save (bool save_at_time_interval = false,
            bool save_at_random_times = false,
            double num_points_to_collect = 1000,
            double burn_in = 5000);

        void init_random_number_generator (std::vector<std::uint_least32_t> random_seeds =
            {582654328, 1065236345, 322147403, 2229968939});

        void set_basal_regulation (std::vector<double> basal_regulation);

        void set_max_freq_regulation (std::vector<double> max_freq_regulation);

        void set_statistics_history_params (unsigned int num_history_statistics = 10000,
            unsigned int burn_in_rxn_count = 10000,
            unsigned int after_burn_in_rxn_count = 10000);

        void set_random_number_generator (RNG &generator);
    };
}
#endif
