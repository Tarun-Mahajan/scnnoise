#include "scnnoise.hpp"
#include "gillespieSSA.hpp"
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <math.h>
// #include "graph.hpp"



namespace ScnnoiseInterface {
  /* function definitions */
  // Constructor
  GillespieSSA::GillespieSSA (int num_genes, std::string gene_filepath,
      std::string molecule_count_filepath,
      std::string count_save_file, bool keep_GRN,
      std::string GRN_filepath, int num_timepoints_save):
    scNNoiSE (num_genes, gene_filepath,
        molecule_count_filepath,
        count_save_file, keep_GRN,
        GRN_filepath, num_timepoints_save) {
  }

  inline double GillespieSSA::sample_time_step (RNG &generator) {
    thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double rand_num = distribution(generator);
    return -log(double(1.0) - rand_num)/total_propensity;
  }

  int GillespieSSA::sample_next_rxn (RNG &generator) {
    thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double rand_num = distribution(generator);
    double selector = total_propensity * (double(1.0) - rand_num);
    int rxn_selected = -1;

    for (std::size_t i = 0; i < rxn_order.size(); ++i) {
      selector -= rxn_order[i].propensity_val;
      if (selector <= 0) {
        rxn_selected = i;
        break;
      }
    }
    gene_type_struct gene_type = gene_type_info[reactions[rxn_order[rxn_selected].gene_id].gene_type];
    if (reactions[rxn_order[rxn_selected].gene_id].molecule_count_cur[gene_type.species_rev_map["mRNA"]] == 0 &&
        rxn_order[rxn_selected].rxn_name == "mRNA decay") {
        std::cout << "error found " <<
        reactions[0].propensity_vals["mRNA decay"] << " " <<
        total_propensity << std::endl;
    }
    sort_reaction(rxn_selected);
    gene_type = gene_type_info[reactions[rxn_order[rxn_selected].gene_id].gene_type];
    if (reactions[rxn_order[rxn_selected].gene_id].molecule_count_cur[gene_type.species_rev_map["mRNA"]] == 0 &&
        rxn_order[rxn_selected].rxn_name == "mRNA decay") {
        std::cout << "error found 1 " <<
        reactions[0].propensity_vals["mRNA decay"] << " " <<
        total_propensity << std::endl;
    }
    // if (rxn_selected > 0) {
    //   std::swap(rxn_order[rxn_selected - 1], rxn_order[rxn_selected]);
    //   rxn_selected -= 1;
    // }
    return rxn_selected;
  }

  // void GillespieSSA::update_fired_reaction_reactants(int gene_selected,
  //   int rxn_index, int &count_not_changed_reactants,
  //   std::vector<bool> &flag_changed_product_count,
  //   std::vector<bool> &GRN_out_changed,
  //   std::vector<int> &rxn_selected_reactants,
  //   std::vector<int> &rxn_selected_products,
  //   std::vector<int> &rxn_selected_reactants_stoichio,
  //   std::vector<int> &rxn_selected_products_stoichio) {
  //     for (auto r = rxn_selected_reactants.begin(); r != rxn_selected_reactants.end(); ++r) {
  //       int reactant_index = std::distance(rxn_selected_reactants.begin(), r);
  //       std::vector<int>::iterator it =
  //       std::find(rxn_selected_products.begin(), rxn_selected_products.end(),
  //       *r);
  //       if (it != rxn_selected_products.end()) {
  //         int product_index = std::distance(rxn_selected_products.begin(), it);
  //         flag_changed_product_count[product_index] = true;
  //         int count_change = rxn_selected_products_stoichio[product_index] -
  //           rxn_selected_reactants_stoichio[reactant_index];
  //         if (count_change) {
  //           reactions[gene_selected].molecule_count_cur[*r] += count_change;
  //           std::vector<int>::iterator it1 =
  //           std::find(reactions[gene_selected].GRN_species_OUT.begin(),
  //           reactions[gene_selected].GRN_species_OUT.end(),
  //           *r);
  //           if (it1 != reactions[gene_selected].GRN_species_OUT.end()) {
  //             int out_index = std::distance(reactions[gene_selected].GRN_species_OUT.begin(), it1);
  //             GRN_out_changed[out_index] = true;
  //           }
  //           // if (*r == reactions[gene_selected].GRN_rxn_OUT) {
  //           //   GRN_out_changed = true;
  //           // }
  //         }else{
  //           count_not_changed_reactants += 1;
  //         }
  //       }else{
  //         reactions[gene_selected].molecule_count_cur[*r] -=
  //           rxn_selected_reactants_stoichio[reactant_index];
  //         std::vector<int>::iterator it1 =
  //         std::find(reactions[gene_selected].GRN_species_OUT.begin(),
  //         reactions[gene_selected].GRN_species_OUT.end(),
  //         *r);
  //         if (it1 != reactions[gene_selected].GRN_species_OUT.end()) {
  //           int out_index = std::distance(reactions[gene_selected].GRN_species_OUT.begin(), it1);
  //           GRN_out_changed[out_index] = true;
  //         }
  //       }
  //     }
  //   }

    std::vector<std::string> GillespieSSA::update_fired_reaction (int rxn_selected) {
        /*
        Find the gene and reaction type for the selected reaction channel. Extract the reactants,
        products and their stoichiometric coefficients for the selected reaction channel.
        */
        int gene_selected = rxn_order[rxn_selected].gene_id;
        std::string rxn_name = rxn_order[rxn_selected].rxn_name;
        std::string gene_type = reactions[gene_selected].gene_type;
        gene_type_struct gene_info = gene_type_info[gene_type];
        int rxn_index = gene_info.rxn_rev_map[rxn_name];
        std::vector<std::string> GRN_species_OUT = reactions[gene_selected].GRN_species_OUT;
        std::vector<std::string> GRN_out_changed;

        // Update fired reaction reactants
        for (auto const &reactants_ : gene_info.rxns[rxn_name].reactants_stoichio) {
            double reactant_factor =
                stoichio_factors[gene_selected].rxns[rxn_name].reactants_factors[reactants_.first];
            int reactant_id = gene_info.species_rev_map[reactants_.first];
            auto it = gene_info.rxns[rxn_name].products_stoichio.find(reactants_.first);
            if (it == gene_info.rxns[rxn_name].products_stoichio.end()) {
                // if (reactant_id == 0 &&
                //     reactions[gene_selected].molecule_count_cur[reactant_id] < int (reactants_.second * reactant_factor)) {
                //     std::cout << "error found " <<
                //     reactions[gene_selected].molecule_count_cur[reactant_id] << " " <<
                //     int (reactants_.second * reactant_factor) << " " <<
                //     reactions[gene_selected].propensity_vals[rxn_name] << " " <<
                //     rxn_order[rxn_selected].propensity_val << " " <<
                //     total_propensity << std::endl;
                // }
                reactions[gene_selected].molecule_count_cur[reactant_id] -=
                    int (reactants_.second * reactant_factor);
                if (keep_GRN) {
                    auto it_out = std::find(GRN_species_OUT.begin(), GRN_species_OUT.end(),
                        reactants_.first);
                    if (it_out != GRN_species_OUT.end()) {
                        GRN_out_changed.push_back(reactants_.first);
                    }
                }
            }else{
                double product_factor =
                    stoichio_factors[gene_selected].rxns[rxn_name].products_factors[reactants_.first];
                int stoichio_diff = int (it->second * product_factor -
                    reactants_.second * reactant_factor);
                if (stoichio_diff != 0) {
                    reactions[gene_selected].molecule_count_cur[reactant_id] +=
                        stoichio_diff;
                    if (keep_GRN) {
                        auto it_out = std::find(GRN_species_OUT.begin(), GRN_species_OUT.end(),
                            reactants_.first);
                        if (it_out != GRN_species_OUT.end()) {
                            GRN_out_changed.push_back(reactants_.first);
                        }
                    }
                }
            }
        }

        // Update fired reaction products
        for (auto const &products_ : gene_info.rxns[rxn_name].products_stoichio) {
            double product_factor =
                stoichio_factors[gene_selected].rxns[rxn_name].products_factors[products_.first];
            int product_id = gene_info.species_rev_map[products_.first];
            auto it = gene_info.rxns[rxn_name].reactants_stoichio.find(products_.first);
            if (it == gene_info.rxns[rxn_name].reactants_stoichio.end()) {
                reactions[gene_selected].molecule_count_cur[product_id] +=
                    int (products_.second * product_factor);
                if (keep_GRN) {
                    auto it_out = std::find(GRN_species_OUT.begin(), GRN_species_OUT.end(),
                        products_.first);
                    if (it_out != GRN_species_OUT.end()) {
                        GRN_out_changed.push_back(products_.first);
                    }
                }
            }
        }
        return GRN_out_changed;
    }

    void GillespieSSA::update_dependent_count_propensity (int rxn_selected,
        std::vector<std::string> &GRN_out_changed) {

        /*
        Update propensity for dependent reactions belonging to the same gene as the
        reaction selected for firing
        */
        int gene_selected = rxn_order[rxn_selected].gene_id;
        std::string rxn_name = rxn_order[rxn_selected].rxn_name;
        std::string gene_type = reactions[gene_selected].gene_type;
        gene_type_struct gene_info = gene_type_info[gene_type];
        int rxn_index = gene_info.rxn_rev_map[rxn_name];
        std::vector<int> rxn_selected_children;
        gene_info.gene_rxn_dependency[0].find_children(rxn_index, rxn_selected_children);
        std::vector<int> rxn_selected_children_GRN;
        if (keep_GRN) {
            network[0].find_children(gene_selected, rxn_selected_children_GRN);
        }
        // rxn_selected_reactants.reserve(reactions[gene_selected].molecule_count_cur.size());
        // std::vector<int> rxn_selected_reactants_stoichio;
        // rxn_selected_reactants_stoichio.reserve(reactions[gene_selected].molecule_count_cur.size());
        for (auto const &rxn : rxn_selected_children) {
            std::string rxn_name_cur = gene_info.rxn_map[rxn];
            total_propensity -= reactions[gene_selected].propensity_vals[rxn_name_cur];
            double new_propensity = compute_propensity (gene_map[gene_selected],
                rxn_name_cur);
            total_propensity += new_propensity;
            reactions[gene_selected].propensity_vals[rxn_name_cur] = new_propensity;
            rxn_order[rxn_order_map[gene_map[gene_selected]][rxn_name_cur]].propensity_val =
                new_propensity;
        }


        /*
        Update propensity for dependent reactions belonging to genes
        downstream of the gene related to the reaction selected for firing
        */
        if (GRN_out_changed.size() != 0 && keep_GRN) {
            unsigned int children_counter = 0;
            for (auto const &dest : rxn_selected_children_GRN) {
                std::string species_OUT =
                    gene_info.species_map[network[0].edge_rxn_params[gene_selected][children_counter].species_OUT];
                std::string gene_type_cur = reactions[dest].gene_type;
                gene_type_struct gene_info_cur = gene_type_info[gene_type_cur];
                std::string rxn_IN =
                    gene_info_cur.rxn_map[network[0].edge_rxn_params[gene_selected][children_counter].rxn_IN];
                auto it = std::find(GRN_out_changed.begin(),
                    GRN_out_changed.end(),
                    species_OUT);
                if (it != GRN_out_changed.end()) {
                    total_propensity -= reactions[dest].propensity_vals[rxn_IN];
                    double new_propensity = compute_propensity (gene_map[dest],
                        rxn_IN);
                    total_propensity += new_propensity;
                    reactions[dest].propensity_vals[rxn_IN] = new_propensity;
                    rxn_order[rxn_order_map[gene_map[dest]][rxn_IN]].propensity_val =
                        new_propensity;
                }
                ++children_counter;
            }
        }
    }

<<<<<<< HEAD
  void GillespieSSA::update_dependent_count_propensity (int rxn_selected,
    std::vector<bool> &GRN_out_changed) {

      /*
      Update propensity for dependent reactions belonging to the same gene as the
      reaction selected for firing
      */
      int gene_selected = rxn_order[rxn_selected].gene_id;
      int rxn_index = rxn_order[rxn_selected].rxn_type;
      int gene_type = reactions[gene_selected].gene_type;
      std::vector<int> rxn_selected_children;
      rxn_selected_children.reserve(gene_rxn_dependency[gene_type].get_size());
      gene_rxn_dependency[gene_type].find_children(rxn_index, rxn_selected_children);
      std::vector<int> rxn_selected_reactants;
      rxn_selected_reactants.reserve(reactions[gene_selected].molecule_count_cur.size());
      std::vector<int> rxn_selected_reactants_stoichio;
      rxn_selected_reactants_stoichio.reserve(reactions[gene_selected].molecule_count_cur.size());
      for (auto &rxn : rxn_selected_children) {
        total_propensity -= reactions[gene_selected].rxns[rxn].propensity_val; //
        rxn_selected_reactants =
          reactions[gene_selected].rxns[rxn].reactants;
        rxn_selected_reactants_stoichio =
          reactions[gene_selected].rxns[rxn].reactants_stoichio;
        double new_propensity = reactions[gene_selected].rxns[rxn].rxn_rate;
        
        for (auto r = rxn_selected_reactants.begin(); r != rxn_selected_reactants.end(); ++r) {
          
          int reactant_index = std::distance(rxn_selected_reactants.begin(), r);
          int reactant_stoichio = rxn_selected_reactants_stoichio[reactant_index];
        
          if (reactions[gene_selected].molecule_count_cur[*r] < reactant_stoichio){
          new_propensity *= 0;
          }else{
          for(int rs = 0; rs<reactant_stoichio; ++rs)
            new_propensity *= (reactions[gene_selected].molecule_count_cur[*r]-rs);
          }
        }
        std::vector<int>::iterator it1 =
        std::find(reactions[gene_selected].GRN_rxn_IN.begin(),
        reactions[gene_selected].GRN_rxn_IN.end(),
        rxn);

        

        
        if (it1 != reactions[gene_selected].GRN_rxn_IN.end()) {
          // int rxn_index =
          // std::distance(reactions[gene_selected].GRN_rxn_IN.begin(), *it1);
          new_propensity *= regulation_function(gene_selected, rxn);
        }

        

        // if (rxn == reactions[gene_selected].GRN_rxn_IN) {
        //   new_propensity *= regulation_function(gene_selected);
        // }
        reactions[gene_selected].rxns[rxn].propensity_val = new_propensity;
        
        for (auto &r : rxn_order) {
          if ((r.gene_id == gene_selected) && (r.rxn_type == rxn)) {
            r.propensity_val = new_propensity;
            break;
          }
        }
        total_propensity += reactions[gene_selected].rxns[rxn].propensity_val;
        
      }
      

      /*
      Update propensity for dependent reactions belonging to genes
      downstream of the gene related to the reaction selected for firing
      */
      for (int g_index = 0; g_index < GRN_out_changed.size(); ++g_index) {
        if (GRN_out_changed[g_index]) {
          std::vector<int> gene_children;
          std::vector<edge_rxn_struct> gene_children_edge_info;
          network[0].find_children(gene_selected, gene_children);
          network[0].find_children_edge_info(gene_selected, gene_children_edge_info);
          // for (auto &g : gene_children) {
          for (int g = 0; g < gene_children.size(); ++g) {

            // int rxn = reactions[gene_selected[g]].GRN_rxn_IN;
            int rxn = gene_children_edge_info[g].rxn_IN;
            int out_species = gene_children_edge_info[g].species_OUT;  
            if (out_species == reactions[gene_selected].GRN_species_OUT[g_index]) {
              total_propensity -= reactions[gene_children[g]].rxns[rxn].propensity_val;
              rxn_selected_reactants =
                reactions[gene_children[g]].rxns[rxn].reactants;
              rxn_selected_reactants_stoichio =
                reactions[gene_children[g]].rxns[rxn].reactants_stoichio;

              double new_propensity = reactions[gene_children[g]].rxns[rxn].rxn_rate;
              
              for (auto r = rxn_selected_reactants.begin(); r != rxn_selected_reactants.end(); ++r) {
                int reactant_index = std::distance(rxn_selected_reactants.begin(), r);
                int reactant_stoichio = rxn_selected_reactants_stoichio[reactant_index];
                if (reactions[gene_selected].molecule_count_cur[*r] < reactant_stoichio){
                  new_propensity *= 0;
                }else{
                for(int rs = 0; rs<reactant_stoichio; ++rs)
                  new_propensity *= (reactions[gene_selected].molecule_count_cur[*r]-rs);
                }
              }
              
              new_propensity *= regulation_function(gene_children[g], rxn);
                           
          
        
              reactions[gene_children[g]].rxns[rxn].propensity_val = new_propensity;
              for (auto &r : rxn_order) {
                if ((r.gene_id == gene_children[g]) && (r.rxn_type == rxn)) {
                  r.propensity_val = new_propensity;
                  break;
                }
              }
 
              total_propensity += reactions[gene_children[g]].rxns[rxn].propensity_val;
=======
    void GillespieSSA::update_molecule_count_history (int &num_history, int &num_save_loop,
        bool simulation_ended) {
        if ((num_history == num_timepoints_save) || simulation_ended) {
            int num_timepoints_save_cur;
            if (simulation_ended) {
                num_timepoints_save_cur = num_history;
            }else{
                num_timepoints_save_cur = num_timepoints_save;
            }
            num_history = 0;
            num_save_loop += 1;
            if (save_timeseries) {
                std::ofstream outfile;
                outfile.open(count_save_file, std::ios_base::app);
                for (int id_time = 0; id_time < num_timepoints_save_cur; ++id_time) {
                    outfile << time_history[(num_save_loop - 1)*num_timepoints_save + id_time] << ',';
                    for (int gene = 0; gene < num_genes; ++gene) {
                        for (int species = 0; species < reactions[gene].molecule_count_cur.size(); ++species) {
                            if (species == reactions[gene].molecule_count_cur.size() - 1 && gene == num_genes - 1) {
                                outfile << molecule_count_history[gene][species][id_time];
                            }else{
                                outfile << molecule_count_history[gene][species][id_time] << ',';
                            }
                        }
                    }
                    outfile << '\n';
                }
                outfile.close();
            }else{
                if (simulation_ended) {
                    std::ofstream outfile;
                    outfile.open(count_save_file, std::ios_base::app);
                    outfile << time_history.back() << ",";
                    for (int gene = 0; gene < num_genes; ++gene) {
                        for (int species = 0; species < reactions[gene].molecule_count_cur.size(); ++species) {
                            if (species == reactions[gene].molecule_count_cur.size() - 1 && gene == num_genes - 1) {
                                outfile << reactions[gene].molecule_count_cur[species];
                            }else{
                                outfile << reactions[gene].molecule_count_cur[species] << ',';
                            }
                        }
                    }
                    outfile << '\n';
                    outfile.close();
                }
            }
        }else{

        }

        for (int gene = 0; gene < num_genes; ++gene) {
          for (int species = 0; species < reactions[gene].molecule_count_cur.size(); ++species) {
            molecule_count_history[gene][species][num_history] =
              reactions[gene].molecule_count_cur[species];
          }
        }
        ++num_history;
    }

    void GillespieSSA::start_molecule_count_history_file () {
        std::ofstream outfile;
        outfile.open(count_save_file);
        outfile << "time" << ",";
        for (int gene = 0; gene < num_genes; ++gene) {
          for (int species = 0; species < reactions[gene].molecule_count_cur.size(); ++species) {
              std::string gene_type = reactions[gene].gene_type;
              std::string connector = ":";
              std::string gene_species_name = gene_map[gene] + connector +
                gene_type_info[gene_type].species_map[species];
            if (species == reactions[gene].molecule_count_cur.size() - 1 && gene == num_genes - 1) {
                outfile << gene_species_name;
            }else{
                outfile << gene_species_name << ",";
>>>>>>> origin/develop
            }
          }
        }
        outfile << '\n';
        outfile << time_history.back() << ",";
        for (int gene = 0; gene < num_genes; ++gene) {
          for (int species = 0; species < reactions[gene].molecule_count_cur.size(); ++species) {
              std::string gene_type = reactions[gene].gene_type;
              std::string connector = ":";
              std::string gene_species_name = gene_map[gene] + connector +
                gene_type_info[gene_type].species_map[species];
            if (species == reactions[gene].molecule_count_cur.size() - 1 && gene == num_genes - 1) {
                outfile << reactions[gene].molecule_count_cur[species];
            }else{
                outfile << reactions[gene].molecule_count_cur[species] << ",";
            }
          }
        }
        outfile << '\n';
        outfile.close();
<<<<<<< HEAD
      }
    }
    for (int gene=0; gene < num_genes; ++gene) {
      for (int species=0; species < num_species_gene_type[gene]; ++species) {
        molecule_count_history[gene][species][num_history] =
          reactions[gene].molecule_count_cur[species];
      }
    }
  }

  void GillespieSSA::simulate (int sim_time) {
    total_propensity = 0;
    compute_total_propensity();
    std::random_device rd;
    std::vector<std::uint_least32_t> rd_seeds = {rd(), rd(), rd(), rd()};
    std::seed_seq sd(rd_seeds.begin(), rd_seeds.end());
    thread_local static RNG generator{sd};
    bool stop_sim = false;
    int num_history = 0;
    int num_save_loop = 0;

    while (!stop_sim) {
      
      num_history += 1;
      double next_time_step = sample_time_step(generator);
      
      int next_rxn = sample_next_rxn(generator);
      
      std::vector<bool> GRN_out_changed;
      double total_time = std::accumulate(time_history.begin(), time_history.end(),
        decltype(time_history)::value_type(0)) + next_time_step;
      
      if (total_time < sim_time) {
          update_cell_cycle_state(total_time);
          
          GRN_out_changed = update_fired_Reaction(next_rxn);
          update_dependent_count_propensity(next_rxn, GRN_out_changed);
          time_history.push_back(next_time_step);
          update_molecule_count_history(num_history, num_save_loop);
        
      }else{
        stop_sim = true;
        num_history = 0;
        num_save_loop += 1;
        if (save_timeseries) {
          std::ofstream outfile;
          outfile.open(count_save_file, std::ios_base::app);
          for (int id_time = 0; id_time < time_history.size()-(num_save_loop - 1)*num_timepoints_save; ++id_time) {
            outfile << time_history[(num_save_loop - 1)*num_timepoints_save + id_time] << ',';
            for (int gene = 0; gene < num_genes; ++gene) {
              for (int species = 0; species < num_species_gene_type[gene]; ++species) {
                outfile << molecule_count_history[gene][species][id_time] << ',';
              }
            }
            outfile << '\n';
          }
          outfile.close();
        }
        time_history.clear();
        time_history.push_back(0);
        std::cout<<time_history.size()<<std::endl;
        /**for (int gene = 0; gene < num_genes; ++gene) {
          for (int species = 0; species < num_species_gene_type[gene]; ++species) {
            molecule_count_history[gene][species].clear();
            molecule_count_history[gene][species].push_back(reactions[gene].molecule_count_cur[species]);
          }
        }*/

      }
=======
    }

    void GillespieSSA::simulate (std::vector<unsigned int> random_seeds) {
        time_history.clear();
        time_history.push_back(0);
        compute_total_propensity();
        start_molecule_count_history_file();

        std::random_device rd;
        // std::vector<std::uint_least32_t> rd_seeds = {rd(), rd(), rd(), rd()};
        // std::vector<std::uint_least32_t> rd_seeds =
        //     {random_seeds[0], random_seeds[1], random_seeds[2], random_seeds[3]};
        std::vector<std::uint_least32_t> rd_seeds(random_seeds.size());
        for (std::size_t b = 0; b < random_seeds.size(); ++b) {
            rd_seeds[b] = random_seeds[b];
        }
        std::seed_seq sd(rd_seeds.begin(), rd_seeds.end());
        thread_local static RNG generator{sd};
        bool stop_sim = false;
        bool simulation_ended = false;
        int num_history = 0;
        int num_save_loop = 0;
        double total_time = 0;
        double cur_time = 0;
        update_molecule_count_history(num_history, num_save_loop,
            simulation_ended);
        init_cell_cycle_state (generator, cur_time);

        while (!stop_sim) {
            // if (reactions[0].molecule_count_cur[0] == 0 && reactions[0].propensity_vals["mRNA decay"] > 0) {
            //     std::cout << "error found0time " <<
            //     reactions[0].propensity_vals["mRNA decay"] << " " <<
            //     total_propensity << " time = " << total_time << std::endl;
            // }
            double next_time_step = sample_time_step(generator);
            // if (reactions[0].molecule_count_cur[0] == 0 && reactions[0].propensity_vals["mRNA decay"] > 0) {
            //     std::cout << "error found0time1 " <<
            //     reactions[0].propensity_vals["mRNA decay"] << " " <<
            //     total_propensity << std::endl;
            // }
            std::vector<std::string> GRN_out_changed;
            cur_time = total_time;
            total_time += next_time_step;
            // std::cout << "reached here 1 = " << std::endl;
            if (total_time < max_time) {
                // if (reactions[0].molecule_count_cur[0] == 0 && reactions[0].propensity_vals["mRNA decay"] > 0) {
                //     std::cout << "error found0 " <<
                //     reactions[0].propensity_vals["mRNA decay"] << " " <<
                //     total_propensity << std::endl;
                // }
                update_cell_cycle_state(total_time, cur_time, generator);
                // if (reactions[0].molecule_count_cur[0] == 0 && reactions[0].propensity_vals["mRNA decay"] > 0) {
                //     std::cout << "error found01 " <<
                //     reactions[0].propensity_vals["mRNA decay"] << " " <<
                //     total_propensity << std::endl;
                // }
                int next_rxn = sample_next_rxn(generator);
                // std::cout << "reached here 2 = " << std::endl;
                GRN_out_changed = update_fired_reaction(next_rxn);
                update_dependent_count_propensity(next_rxn, GRN_out_changed);
                // if (reactions[0].molecule_count_cur[0] == 0 && reactions[0].propensity_vals["mRNA decay"] > 0) {
                //     std::cout << "error found3 " <<
                //     reactions[0].propensity_vals["mRNA decay"] << " " <<
                //     total_propensity << std::endl;
                // }
                time_history.push_back(next_time_step);
                update_molecule_count_history(num_history, num_save_loop,
                    simulation_ended);
                // if (reactions[0].molecule_count_cur[0] == 0 && reactions[0].propensity_vals["mRNA decay"] > 0) {
                //     std::cout << "error foundlast " <<
                //     reactions[0].propensity_vals["mRNA decay"] << " " <<
                //     total_propensity << std::endl;
                // }
            }else{
                simulation_ended = true;
                update_molecule_count_history(num_history, num_save_loop,
                    simulation_ended);
                stop_sim = true;
            }
        }
>>>>>>> origin/develop
    }
}
