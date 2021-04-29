// scNNoiSE header file
#ifndef SCNNOISE
#define SCNNOISE
class scNNoiSE {
protected:
  /* data */
  /********************************************//**
  \brief Number of chemical reaction channels.

  This is the total number of unique chemical
  raction channels that are present in the  user
  provided chemical reaction network.
   ***********************************************/
  int num_rxns;

  /********************************************//**
  \brief Number of molecular species.

  This is the total number of unique molecular
  species that are present in the user provided
  chemical reaction network.
   ***********************************************/
  int num_species;

  /********************************************//**
   \brief Vector of current molecule counts for all species.

   This stores the current molecule counts for all
   species in the chemical reaction network along
   the simulation trajectory.
   ***********************************************/
  std::vector<int> molecule_count_cur;

   /********************************************//**
    \brief Reaction search order

    Order in which reaction propensities are added
    for reaction channel selection at each step of the
    simulation.
    ***********************************************/
  std::vector<int> rxn_order;

  /********************************************//**
   \brief Gene regulatory network.

   Gene regulatory network is stored in a vector of
   type GRN, which is a class for gene regulatory
   networks.
   ***********************************************/
  std::vector<GRN> network;

  /********************************************//**
   \brief Vector for node type in the GRN

   A vector to identify whether nodes in the GRN are
   genes or miRNA. Further, it also identifies
   whether the gene is producing only mRNA or both
   nascent and mature mRNA. Each node is labeled by
   a string. The potential labels are as follows:
   'GMP' : gene transcribes mRNA which translates into
   protein
   'GNMP' : gene transcribes nascent mRNA which processes
   into mature mRNA. Mature mRNA is translated into
   protein.
   'Mi' : gene transcribes nascent miRNA, which is processed
   into mature miRNA.
   ***********************************************/
  std::vector<std:string> node_type;

  /********************************************//**
   \brief Number of splicing variants

   Number of different variant mRNA species generated
   by alternative splicing.
   ***********************************************/
  std::vector<int> num_splice_variants;

  /********************************************//**
   \brief Map to store all the reaction channels

   A map of maps that stores information about the reaction
   channels. For each reaction channel, the different
   keys are--'reactants', 'products', 'reactant stoichiometry',
   'product stoichiometry'.
   ***********************************************/
  std::map<int, std::map<std:string, std::vector<int>>> reactions;

  /********************************************//**
   \brief Vector to store reaction rates.
   ***********************************************/
  std::vector<double> rxn_rates;

  // /********************************************//**
  //  \brief Vector for status in GRN
  //
  //  Vector to store status of the molecular species--
  //  whether the species participates in the GRN or not.
  //  ***********************************************/
  // std::vector<bool> in_GRN;

  /********************************************//**
   \brief Vector for order of GRN nodes among species

   The order in which the species occur in the nodes
   in the GRN. If a molecular species is not part of
   the GRN, it is assigned a value of -1.
   ***********************************************/
  std::vector<int> species_order_GRN;

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

public:
  /********************************************//**
   \brief Constructor for scNNoiSE.

   \param[in] num_rxns Number of channels in the
              list of chemical reaction channels.
   \param[in] num_species number of molecular species
              in the system.
   ***********************************************/
  scNNoiSE (int num_rxns, int num_species);

  /********************************************//**
   \brief Destructor for scNNoiSE.
   ***********************************************/
  virtual ~scNNoiSE ();

  /********************************************//**
   \brief Initialize counts for molecular species.

   Function to initialize counts for a single molecular
   species.
   \param[in] molecule_id Integer id for the molecular species
              to initialize.
   \param[in] molecule_count Integer count for the molecular
              species to initialize.
   ***********************************************/
  void init_molecular_count (int molecule_id, int molecule_count);

  /********************************************//**
   \brief Add node-type for a species in the GRN.

   Nodes in the GRN are genes or miRNA. Node-type
   is a label which identifies whether a node is a
   gene or an miRNA. Further, it also identifies
   whether the gene is producing only mRNA or both
   nascent and mature mRNA. Any node is labeled by
   a string. The potential labels are as follows:
   'GMP' : gene transcribes mRNA which translates into
   protein
   'GNMP' : gene transcribes nascent mRNA which processes
   into mature mRNA. Mature mRNA is translated into
   protein.
   'Mi' : gene transcribes nascent miRNA, which is processed
   into mature miRNA.
   \param[in] node_id Integer id for the molecular species
              in the GRN.
   \param[in] str_type String identifying the node-type for the
              molecular species identified by node_id.
   ***********************************************/
  void add_GRNnode_type (int node_id, std::string str_type);

  /********************************************//**
   \brief Add number of alternative splicing variants.

   Function to add number of alternative splicing (AS) variants
   for a molecular species.
   \param[in] gene_id Integer id for the gene for which
              a molecular species is present in the GRN.
   \param[in] num_AS number of AS variants for gene identified by
              node_id.
   ***********************************************/
  void add_num_splice_variants (int gene_id, int num_AS);

  /********************************************//**
   \brief Add reactant stoichiometry for a reaction.

   Function to add reactant stoichiometry for a reaction
   in the chemical system.
   \param[in] rxn_id Integer id for the reaction in the list
              of chemical reaction channels.
   \param[in] rxn_reactants vector of integer molecular species ids
              which are reactants for reaction rxn_id .
   \param[in] reactants_stoichio vector of integer
              stoichoimetry coefficients for reactants in
              rxn_reactants.
   ***********************************************/
  void add_rxn_reactants (int rxn_id, std::vector<int> rxn_reactants,
  std::vector<int> reactants_stoichio);

  /********************************************//**
   \brief Add product stoichiometry for a reaction.

   Function to add product stoichiometry for a reaction
   in the chemical system.
   \param[in] rxn_id Integer id for the reaction in the list
              of chemical reaction channels.
   \param[in] rxn_products vector of integer molecular species ids
              which are products for reaction rxn_id .
   \param[in] products_stoichio vector of integer
              stoichoimetry coefficients for products in
              rxn_products.
   ***********************************************/
  void add_rxn_products (int rxn_id, std::vector<int> rxn_products,
  std::vector<int> products_stoichio);

  /********************************************//**
   \brief Add reaction rate for a reaction in the
          list of chemical reaction channels.

   Function to add reaction rate for a reaction
   in the list of chemical reaction channels.
   \param[in] rxn_id Integer id for the reaction in the list
              of chemical reaction channels.
   \param[in] rxn_rate reaction rate for the reaction
              rxn_id.
   ***********************************************/
  void add_rxn_rates (int rxn_id, double rxn_rate);

  /********************************************//**
   \brief Add species order in GRN.

   Function to add order in which a molecular species
   occurs in the node list of GRN.
   \param[in] molecule_id Integer id for the molecular
              species.
   \param[in] GRN_order of molecule species
              molecule_id in the GRN. If the molecular
              species is not present in the GRN, it is
              assigned a value of -1.
   ***********************************************/
  void add_species_order_GRN (int molecule_id, int GRN_order);

  /********************************************//**
   \brief Simulating stochastic gene expression.

   A pure virtual function for simulating stochastic gene
   expression. Needs to be overridden in any derived class.
   ***********************************************/
  virtual void simulate () = 0;
};
#endif
