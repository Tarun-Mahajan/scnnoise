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



namespace ScnnoiseInterface {
  /********************************************//**
   \brief Struct to store information for each reaction in
          the reaction order list.
   ***********************************************/
    struct rxn_order_struct {
        int gene_id;
        std::string rxn_name;
        double propensity_val;
    };

    struct rxn_struct {
        // std::vector<int> reactants;
        // std::vector<int> products;
        std::map<std::string, int> reactants_stoichio;
        std::map<std::string, int> products_stoichio;
        // double rxn_rate;
        // double propensity_val;
    };

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

    struct stoichio_factor_species_struct {
        // std::vector<int> reactants;
        // std::vector<int> products;
        std::map<std::string, double> reactants_factors;
        std::map<std::string, double> products_factors;
        // double rxn_rate;
        // double propensity_val;
    };

    struct stoichio_factor_struct {
        // std::vector<int> reactants;
        // std::vector<int> products;
        std::map<std::string, stoichio_factor_species_struct> rxns;
        // double rxn_rate;
        // double propensity_val;
    };

  /********************************************//**
   \brief Struct to store information for all reaction
          channels for a gene.

   The data members of the struct are:

   1. gene_type -- An integer to identify whether the node is a
   gene or an miRNA. Further, it also identifies
   whether the gene is producing only mRNA or both
   nascent and mature mRNA. The node is labeled by
   an integer representing a unique string.
   The potential string labels and their corresponding integer labels
   are as follows:
   0 : 'GMP' : gene transcribes mRNA which translates into
   protein. 4 and 6 species and reactions in the chemical reaction
   network, respectively.
   1 : 'GNMP' : gene transcribes nascent mRNA which processes
   into mature mRNA. Mature mRNA is translated into
   protein. For \f$n\f$ alternatively spliced mature mRNA,
   there are \f$2 + 2n\f$ and \f$2 + 4n\f$ chemical species and reactions in the
   chemical reaction network, respectively.
   2 : 'Mi' : gene transcribes nascent miRNA, which is processed
   into mature miRNA.  4 and 5 species and reactions in the chemical
   reaction network, respectively.
   3 : 'GNMP-Mihost' : gene transcribes nascent mRNA which processes
   into mature mRNA and nascent miRNA. Mature mRNA is translated into
   protein. Nascent miRNA is processed into mature miRNA.
   For \f$n\f$ alternatively spliced mature mRNA,
   there are \f$4 + 2n\f$ and \f$4 + 4n\f$ chemical species and reactions in the
   chemical reaction network, respectively.

   2. num_splice_variants -- number of alternatively spliced (AS) mature mRNA
   species.

   3. GRN_rxn_IN -- vector of integer IDs for chemical reactions for which reaction
   rate is affected by the GRN GRN_species_OUT

   4. GRN_species_OUT -- vector of integer IDs for chemical species which cause
   regulation in the GRN

   5. rxn_type -- a vector to identify reaction type for all the reaction
   channels for the gene. Each integer in the vector
   corresponds to a single reaction channel.
   For 'GMP'--
   0 : promoter activation-- promoter switches from the off to the on
                           state
   1 : promoter inactivation -- promoter switches from the on to the off
                              state
   2 : transcription of mRNA
   3 : mRNA degradation
   4 : protein translation from mRNA
   5 : protein degradation

   In the following, \f$n\f$ is the number of alternatively spliced (AS)
   mRNAs
   For 'GNMP' --
   0 : promoter activation-- promoter switches from the off to the on
       state
   1 : promoter inactivation -- promoter switches from the on to the off
       state
   2 : transcription of nascent mRNA
   3 : maturation of nascent mRNA into ASed mature mRNA 1
   ...
   2 + n : maturation of nascent mRNA into ASed mature mRNA n
   3 + n : degradation ASed mature mRNA 1
   ...
   2 + 2n : degradation ASed mature mRNA n
   3 + 2n : protein translation from ASed mature mRNA 1
   ...
   2 + 3n : protein translation from ASed mature mRNA n
   3 + 3n : protein degradation from ASed protein 1
   ...
   2 + 4n : protein degradation from ASed protein n

   For 'Mi' --
   0 : promoter activation-- promoter switches from the off to the on
       state
   1 : promoter inactivation -- promoter switches from the on to the off
       state
   2 : transcription of nascent miRNA
   3 : maturation of nascent miRNA to mature miRNA
   4 : mature miRNA degradation

   In the following, \f$n\f$ is the number of alternatively spliced (AS)
   mRNAs
   For 'GNMP-Mihost' --
   0 : promoter activation-- promoter switches from the off to the on
       state
   1 : promoter inactivation -- promoter switches from the on to the off
       state
   2 : transcription of nascent mRNA
   3 : maturation of nascent mRNA into ASed mature mRNA 1 and nascent
       miRNA-host
   ...
   2 + n : maturation of nascent mRNA into ASed mature mRNA n  and
           nascent miRNA-host
   3 + n : degradation ASed mature mRNA 1
   ...
   2 + 2n : degradation ASed mature mRNA n
   3 + 2n : protein translation from ASed mature mRNA 1
   ...
   2 + 3n : protein translation from ASed mature mRNA n
   3 + 3n : protein degradation from ASed protein 1
   ...
   2 + 4n : protein degradation from ASed protein n
   3 + 4n : maturation of nascent miRNA-host into mature miRNA-host
   4 + 4n : degradation of mature miRNA-host

   6. rxns a vector of structs which store stoichiometry and rate
   information for all the reaction channels for the gene. Each struct
   in the vector corresponds to a single reaction channel.

   7. molecule_count_cur vector of current molecule count for chemical
   species related to the gene. The order of the chemical species in
   vector, depedning on the gene type, is as follows:
   For 'GMP'--
   0 : promoter active state
   1 : promoter inactive state
   2 : mRNA
   3 : protein

   In the following, \f$n\f$ is the number of alternatively spliced (AS)
   mRNAs
   For 'GNMP' --
   0 : promoter active state
   1 : promoter inactive state
   2 : nascent mRNA
   3 : ASed mature mRNA 1
   ...
   2 + n : ASed mature mRNA n
   3 + n : protein from ASed mature mRNA 1
   ...
   2 + 2n : protein from ASed mature mRNA n

   For 'Mi' --
   0 : promoter active state
   1 : promoter inactive state
   2 : nascent miRNA
   3 : mature miRNA

   In the following, \f$n\f$ is the number of alternatively spliced (AS)
   mRNAs
   For 'GNMP-Mihost' --
   0 : promoter active state
   1 : promoter inactive state
   2 : nascent mRNA
   3 : ASed mature mRNA 1
   ...
   2 + n : ASed mature mRNA n
   3 + n : protein from ASed mature mRNA 1
   ...
   2 + 2n : protein from ASed mature mRNA n
   3 + 2n : nascent miRNA-host
   4 + 2n : mature miRNA-host

   ***********************************************/
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

    struct molecule_history_struct {
        std::vector<std::vector<int>> molecule_count;
    };

    class scNNoiSE {
    public:
        typedef std::mt19937 RNG;
    /* data */
        std::vector<RNG> generator;
        /********************************************//**
        \brief Number of chemical reaction channels.

        This is the total number of unique chemical
        reaction channels that are present in the  user
        provided chemical reaction network.
         ***********************************************/
        int num_rxns;

        /********************************************//**
        \brief Number of nodes in GRN.

        This is the total number of unique chemical
        species that are present in the user
        provided GRN.
         ***********************************************/
        int num_genes;

        std::map<std::string, gene_type_struct> gene_type_info;

        std::map<int, std::string> gene_map;

        std::map<std::string, int> gene_rev_map;

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

        /********************************************//**
        \brief Map to store all the reaction channels

        A map of maps that stores information about the reaction
        channels. For each reaction channel, the different
        keys are--'Gene', 'reaction', 'reactants', 'products', 'reactant stoichiometry',
        'product stoichiometry', 'dependent'. For the 'reaction' key, the different reaction types
        are encoded with the following types along with their integer keys:

        For 'GMP'--
        0 : promoter activation-- promoter switches from the off to the on
                             state
        1 : promoter inactivation -- promoter switches from the on to the off
                                state
        2 : transcription of mRNA
        3 : mRNA degradation
        4 : protein translation from mRNA
        5 : protein degradation

        In the following, \f$n\f$ is the number of alternatively spliced (AS) mRNAs
        For 'GNMP' --
        0 : promoter activation-- promoter switches from the off to the on
                             state
        1 : promoter inactivation -- promoter switches from the on to the off
                                state
        2 : transcription of nascent mRNA
        3 : maturation of nascent mRNA into ASed mature mRNA 1
        ...
        2 + n : maturation of nascent mRNA into ASed mature mRNA n
        3 + n : degradation ASed mature mRNA 1
        ...
        2 + 2n : degradation ASed mature mRNA n
        3 + 2n : protein translation from ASed mature mRNA 1
        ...
        2 + 3n : protein translation from ASed mature mRNA n
        3 + 3n : protein degradation from ASed protein 1
        ...
        2 + 4n : protein degradation from ASed protein n

        For 'Mi' --
        0 : promoter activation-- promoter switches from the off to the on
                             state
        1 : promoter inactivation -- promoter switches from the on to the off
                                state
        2 : transcription of nascent miRNA
        3 : maturation of nascent miRNA to mature miRNA
        4 : mature miRNA degradation

        In the following, \f$n\f$ is the number of alternatively spliced (AS) mRNAs
        For 'GNMP-Mihost' --
        0 : promoter activation-- promoter switches from the off to the on
                             state
        1 : promoter inactivation -- promoter switches from the on to the off
                                state
        2 : transcription of nascent mRNA
        3 : maturation of nascent mRNA into ASed mature mRNA 1 and nascent miRNA-host
        ...
        2 + n : maturation of nascent mRNA into ASed mature mRNA n  and nascent miRNA-host
        3 + n : degradation ASed mature mRNA 1
        ...
        2 + 2n : degradation ASed mature mRNA n
        3 + 2n : protein translation from ASed mature mRNA 1
        ...
        2 + 3n : protein translation from ASed mature mRNA n
        3 + 3n : protein degradation from ASed protein 1
        ...
        2 + 4n : protein degradation from ASed protein n
        3 + 4n : maturation of nascent miRNA-host into mature miRNA-host
        4 + 4n : degradation of mature miRNA-host
        ***********************************************/
        // std::map<int, std::map<std:string, std::vector<int>>> reactions;
        // std::vector<std::map<int, std::map<std:string, std::vector<int>>>> reactions;
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

        gene_type_struct create_two_state_reduced_mRNA_type ();

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
            std::string statistics_file = "dummy", bool verbose = true) = 0;

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

        void set_count_rxns_fired (bool count_rxns,
            unsigned int stop_rxn_count = pow(10, 6));

        void update_rxn_count (int rxn_selected, bool &stop_sim, bool &reached_rxn_count);

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
