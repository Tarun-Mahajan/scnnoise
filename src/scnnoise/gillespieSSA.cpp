#include "scnnoise.hpp"
#include "gillespieSSA.hpp"
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cctype>
#include <sstream>
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
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double rand_num = distribution(generator);
    return -log(double(1.0) - rand_num)/total_propensity;
  }

  int GillespieSSA::sample_next_rxn (RNG &generator) {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
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
    if (rxn_selected == -1) {
        std::uniform_int_distribution<int> distribution(0, rxn_order.size() - 1);
        rxn_selected = distribution(generator);
    }
    // gene_type_struct gene_type = gene_type_info[reactions[rxn_order[rxn_selected].gene_id].gene_type];
    // if (reactions[rxn_order[rxn_selected].gene_id].molecule_count_cur[0] == 0 &&
    //     rxn_order[rxn_selected].rxn_name == "mRNA decay") {
    //     std::cout << "error found " <<
    //     reactions[rxn_order[rxn_selected].gene_id].propensity_vals["mRNA decay"] << " " <<
    //     total_propensity << " selector " << selector_tmp << std::endl;
    // }
    sort_reaction(rxn_selected);
    // gene_type = gene_type_info[reactions[rxn_order[rxn_selected].gene_id].gene_type];
    // if (reactions[rxn_order[rxn_selected].gene_id].molecule_count_cur[0] == 0 &&
    //     rxn_order[rxn_selected].rxn_name == "mRNA decay") {
    //     std::cout << "error found 1 " <<
    //     reactions[0].propensity_vals["mRNA decay"] << " " <<
    //     total_propensity << std::endl;
    // }
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

    void GillespieSSA::save_molecule_count_at_interval (double time_prev, double time_next,
        double &points_collected_prev, std::string next_rxn_name) {
        if (time_next > burn_in) {
            if (((time_next - burn_in) >= (points_collected_prev + 1) * time_interval_to_save) &&
                ((time_prev - burn_in) < (points_collected_prev + 1) * time_interval_to_save)) {
                std::ofstream outfile;
                outfile.open(count_save_file, std::ios_base::app);
                outfile << next_rxn_name << ",";
                outfile << get_cur_cell_cycle_state() << ",";
                outfile << time_next << ',';
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
                points_collected_prev += 1;
            }
        }
    }

    void GillespieSSA::save_molecule_count_at_random_times (double time_prev, double time_next,
        unsigned int &which_random_time_saved, int next_rxn, 
        double time_since_last_division) {
        while (((time_next) > random_times_to_save[which_random_time_saved]) &&
            ((time_prev) > random_times_to_save[which_random_time_saved]) &&
            which_random_time_saved < num_points_to_collect - 1) {
                ++which_random_time_saved;
        }
        if (time_next > burn_in && time_next < max_time) {
            if (((time_next) >= random_times_to_save[which_random_time_saved]) &&
                ((time_prev) < random_times_to_save[which_random_time_saved])) {
                std::string cur_cell_cycle_phase = get_cur_cell_cycle_state();
                std::ofstream outfile;
                outfile.open(count_save_file, std::ios_base::app);
                outfile << rxn_order[next_rxn].rxn_name << ",";
                outfile << cur_cell_cycle_phase << ',';
                outfile << time_next << ',';
                outfile << time_since_last_division << ",";
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
                which_random_time_saved += 1;
            }
        }
    }

    void GillespieSSA::update_molecule_count_history (int &num_history, int &num_save_loop,
        bool simulation_ended, double cur_time, std::string cur_rxn, 
        double time_since_last_division) {
        if ((num_history == num_timepoints_save && (save_timeseries || save_timeseries_all)) || simulation_ended) {
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
                if (save_timeseries_all) {
                    outfile.open(count_save_file, std::ios_base::app);
                }else{
                    if (!simulation_ended) {
                        outfile.open(count_save_file);
                        outfile << "reaction" << ",";
                        outfile << "phase" << ",";
                        outfile << "time" << ",";
                        outfile << "time_since_last_division" << ",";
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
                            }
                          }
                        }
                        outfile << '\n';
                    }else{
                        outfile.open(count_save_file, std::ios_base::app);
                    }
                }
                for (int id_time = 0; id_time < num_timepoints_save_cur; ++id_time) {
                    outfile << rxn_history[id_time] << ",";
                    outfile << cell_cycle_phase_history[id_time] << ",";
                    // outfile << time_history[(num_save_loop - 1)*num_timepoints_save + id_time] << ',';
                    outfile << time_history[id_time] << ',';
                    outfile << time_since_last_division << ',';
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
                    outfile << cur_rxn << ",";
                    outfile << get_cur_cell_cycle_state() << ",";
                    outfile << cur_time << ",";
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

        if (save_timeseries || save_timeseries_all) {
            rxn_history[num_history] = cur_rxn;
            time_history[num_history] = cur_time;
            cell_cycle_phase_history[num_history] = get_cur_cell_cycle_state();
            for (int gene = 0; gene < num_genes; ++gene) {
              for (int species = 0; species < reactions[gene].molecule_count_cur.size(); ++species) {
                molecule_count_history[gene][species][num_history] =
                  reactions[gene].molecule_count_cur[species];
              }
            }
            ++num_history;
        }
    }

    void GillespieSSA::start_molecule_count_history_file () {
        std::ofstream outfile;
        outfile.open(count_save_file);
        outfile << "reaction" << ",";
        outfile << "phase" << ",";
        outfile << "time" << ",";
        outfile << "time_since_last_division" << ",";
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
            }
          }
        }
        outfile << '\n';
        // outfile << time_history.back() << ",";
        // for (int gene = 0; gene < num_genes; ++gene) {
        //   for (int species = 0; species < reactions[gene].molecule_count_cur.size(); ++species) {
        //       std::string gene_type = reactions[gene].gene_type;
        //       std::string connector = ":";
        //       std::string gene_species_name = gene_map[gene] + connector +
        //         gene_type_info[gene_type].species_map[species];
        //     if (species == reactions[gene].molecule_count_cur.size() - 1 && gene == num_genes - 1) {
        //         outfile << reactions[gene].molecule_count_cur[species];
        //     }else{
        //         outfile << reactions[gene].molecule_count_cur[species] << ",";
        //     }
        //   }
        // }
        // outfile << '\n';
        outfile.close();
    }

    void GillespieSSA::start_statistics_file (std::string statistics_file) {
        std::ofstream outfile;
        outfile.open(statistics_file);
        std::string connector = ":";
        std::string connector_2 = "_";
        for (int gene = 0; gene < num_genes; ++gene) {
            std::string gene_type = reactions[gene].gene_type;
            for (int species = 0; species < reactions[gene].molecule_count_cur.size(); ++species) {
                std::string gene_species_name = gene_map[gene] + connector +
                gene_type_info[gene_type].species_map[species];

                for (unsigned int gene_2 = 0; gene_2 < num_genes; ++gene_2) {
                    std::string gene_type_2 = reactions[gene_2].gene_type;
                    for (int species_2 = 0; species_2 < reactions[gene_2].molecule_count_cur.size(); ++species_2) {
                        std::string gene_species_name_2 = gene_map[gene_2] + connector +
                        gene_type_info[gene_type_2].species_map[species_2];
                        outfile << "cov_" + gene_species_name + connector_2 +
                            gene_species_name_2 << ",";
                    }
                }

                if (species == reactions[gene].molecule_count_cur.size() - 1 && gene == num_genes - 1) {
                    outfile << "mean_" + gene_species_name << "," <<
                        "variance_" + gene_species_name << "\n";
                }else{
                    outfile << "mean_" + gene_species_name << "," <<
                        "variance_" + gene_species_name << ",";
                }
            }
        }
        outfile.close();
    }

    void GillespieSSA::write_statistics_to_file (std::string statistics_file) {
        std::ofstream outfile;
        outfile.open(statistics_file, std::ios_base::app);
        unsigned int count_ = 0;
        for (int gene_ = 0; gene_ < num_genes; ++gene_) {
            for (int species_ = 0; species_ < reactions[gene_].molecule_count_cur.size(); ++species_) {

                for (unsigned int gene_2 = 0; gene_2 < num_genes; ++gene_2) {
                    for (int species_2 = 0; species_2 < reactions[gene_2].molecule_count_cur.size(); ++species_2) {
                        outfile << running_cov[count_] << ",";
                        count_ += 1;
                    }
                }

                if (species_ == reactions[gene_].molecule_count_cur.size() - 1 && gene_ == num_genes - 1) {
                    outfile << running_mean[gene_][species_] << "," <<
                        running_var[gene_][species_] << "\n";
                }else{
                    outfile << running_mean[gene_][species_] << "," <<
                        running_var[gene_][species_] << ",";
                }
            }
        }
        outfile.close();
    }

    void GillespieSSA::set_size_statistics_containers () {
        running_mean.resize(num_genes);
        running_var.resize(num_genes);

        unsigned int count_ = 0;
        for (unsigned int gene_ = 0; gene_ < num_genes; ++gene_) {
            for (int species_ = 0; species_ < reactions[gene_].molecule_count_cur.size(); ++species_) {
                for (unsigned int gene_2 = 0; gene_2 < num_genes; ++gene_2) {
                    for (int species_2 = 0; species_2 < reactions[gene_2].molecule_count_cur.size(); ++species_2) {
                        count_ += 1;
                    }
                }
            }
            running_mean[gene_].resize(reactions[gene_].molecule_count_cur.size(), 0);
            running_var[gene_].resize(reactions[gene_].molecule_count_cur.size(), 0);
        }

        running_cov.resize(count_, 0);
    }

    void GillespieSSA::set_size_statistics_history_containers () {
        history_statistics.resize(num_genes);
        for (unsigned int gene_ = 0; gene_ < num_genes; ++gene_) {
            history_statistics[gene_].resize(reactions[gene_].molecule_count_cur.size());
            unsigned int mol_size = reactions[gene_].molecule_count_cur.size();
            for (unsigned int mol_ = 0; mol_ < mol_size; ++mol_) {
                history_statistics[gene_][mol_].resize(num_history_statistics, 0);
            }
        }
    }

    void GillespieSSA::upate_running_statistics (double total_time_prev, double next_time_step) {
        double total_time_next = total_time_prev + next_time_step;
        double time_ratio = 0;
        if (total_time_prev > 0) {
            time_ratio = pow(10, log10(total_time_prev) -
                log10(total_time_next));
        }
        unsigned int  count_ = 0;
        for (unsigned int gene_ = 0; gene_ < num_genes; ++gene_) {
            for (int species_ = 0; species_ < reactions[gene_].molecule_count_cur.size(); ++species_) {
                double prev_count = reactions[gene_].molecule_count_cur[species_];
                double prev_mean = running_mean[gene_][species_];
                double next_mean =
                    prev_mean * (time_ratio);
                next_mean +=
                    prev_count * (next_time_step / total_time_next);
                for (unsigned int gene_2 = 0; gene_2 < num_genes; ++gene_2) {
                    for (int species_2 = 0; species_2 < reactions[gene_2].molecule_count_cur.size(); ++species_2) {
                        double prev_count_2 =
                            reactions[gene_2].molecule_count_cur[species_2];
                        double prev_mean_2 = running_mean[gene_2][species_2];
                        double next_mean_2 =
                            prev_mean_2 * (time_ratio);
                        next_mean_2 +=
                            prev_count_2 * (next_time_step / total_time_next);
                        running_cov[count_] += prev_mean * prev_mean_2;
                        running_cov[count_] *=
                            (total_time_prev / total_time_next);
                        running_cov[count_] += (prev_count * prev_count_2) *
                            (next_time_step / total_time_next);
                        running_cov[count_] -= next_mean * next_mean_2;
                        count_ += 1;
                    }
                }

                running_mean[gene_][species_] = next_mean;
                // running_mean[gene_][species_] =
                //     prev_mean * (total_time_prev / total_time_next);
                // running_mean[gene_][species_] +=
                //     prev_count * (next_time_step / total_time_next);
                running_var[gene_][species_] =
                    running_var[gene_][species_] * (time_ratio);
                running_var[gene_][species_] +=
                    (pow(prev_count, 2) * (next_time_step / total_time_next));
                running_var[gene_][species_] +=
                    (pow(prev_mean, 2) * (time_ratio));
                running_var[gene_][species_] -=
                    pow(running_mean[gene_][species_], 2);
            }
        }
    }

    void GillespieSSA::upate_running_statistics_without_time_step (
        double num_rxns_cur) {
        double num_rxns_next = num_rxns_cur + 1;
        double num_rxn_ratio = 0;
        if (num_rxns_cur != 0) {
            num_rxn_ratio = pow(10, log10(num_rxns_cur) -
                log10(num_rxns_next));
        }

        unsigned int  count_ = 0;
        for (unsigned int gene_ = 0; gene_ < num_genes; ++gene_) {
            for (int species_ = 0; species_ < reactions[gene_].molecule_count_cur.size(); ++species_) {
                double prev_count = reactions[gene_].molecule_count_cur[species_];
                double prev_mean = running_mean[gene_][species_];
                double next_mean =
                    prev_mean * num_rxn_ratio;
                next_mean +=
                    prev_count * (1.0 / num_rxns_next);
                for (unsigned int gene_2 = 0; gene_2 < num_genes; ++gene_2) {
                    for (int species_2 = 0; species_2 < reactions[gene_2].molecule_count_cur.size(); ++species_2) {
                        double prev_count_2 =
                            reactions[gene_2].molecule_count_cur[species_2];
                        double prev_mean_2 = running_mean[gene_2][species_2];
                        double next_mean_2 =
                            prev_mean_2 * num_rxn_ratio;
                        next_mean_2 +=
                            prev_count_2 * (1.0 / num_rxns_next);
                        running_cov[count_] += prev_mean * prev_mean_2;
                        running_cov[count_] *=
                            num_rxn_ratio;
                        running_cov[count_] += (prev_count * prev_count_2) *
                            (1.0 / num_rxns_next);
                        running_cov[count_] -= next_mean * next_mean_2;
                        count_ += 1;
                    }
                }

                running_mean[gene_][species_] = next_mean;
                // running_mean[gene_][species_] =
                //     prev_mean * (total_time_prev / total_time_next);
                // running_mean[gene_][species_] +=
                //     prev_count * (next_time_step / total_time_next);
                running_var[gene_][species_] =
                    running_var[gene_][species_] * num_rxn_ratio;
                running_var[gene_][species_] +=
                    (pow(prev_count, 2) * (1.0 / num_rxns_next));
                running_var[gene_][species_] +=
                    (pow(prev_mean, 2) * num_rxn_ratio);
                running_var[gene_][species_] -=
                    pow(running_mean[gene_][species_], 2);
            }
        }
    }

    void GillespieSSA::set_running_probs_containers (unsigned int running_probs_buffer_size,
        std::string filepath_probs, double total_time_norm) {
        this->running_probs_buffer_size = running_probs_buffer_size;
        this->filepath_probs = filepath_probs;
        this->total_time_norm = total_time_norm;
        running_marginal_probs.resize(num_genes);

        unsigned int count_ = 0;
        for (unsigned int gene_ = 0; gene_ < num_genes; ++gene_) {
            for (int species_ = 0; species_ < reactions[gene_].molecule_count_cur.size(); ++species_) {
                for (unsigned int gene_2 = 0; gene_2 < num_genes; ++gene_2) {
                    for (int species_2 = 0; species_2 < reactions[gene_2].molecule_count_cur.size(); ++species_2) {
                        if ((gene_ == gene_2) && (species_ == species_2)) {

                        } else {
                            count_ += 1;
                        }
                    }
                }
            }
            running_marginal_probs[gene_].resize(reactions[gene_].molecule_count_cur.size());
        }

        running_joint_probs.resize(count_);
    }

    void GillespieSSA::output_probs_to_file (bool if_first_write) {
        std::string connector_ = "_";
        std::string file_ext = ".csv";
        unsigned int  count_ = 0;
        for (unsigned int gene_ = 0; gene_ < num_genes; ++gene_) {
            for (int species_ = 0; species_ < reactions[gene_].molecule_count_cur.size(); ++species_) {
                std::string filepath = filepath_probs + connector_ +
                    std::to_string(gene_) + connector_ +
                    std::to_string(species_) + file_ext;

                if (!if_first_write) {
                    std::ifstream prob_state_file(filepath);
                    std::string row_text;
                    unsigned int count_key;
                    double prob_;
                    std::string word;
                    unsigned int line_counter = 0;
                    while (std::getline(prob_state_file, row_text)) {
                        if (line_counter > 0) {
                            std::istringstream str_stream(row_text);

                            unsigned int id_counter = 0;
                            while (std::getline(str_stream, word, ',')) {
                                switch(id_counter) {
                                    case 0:
                                        {
                                            count_key = std::stoi(word);
                                            break;
                                        }
                                    case 1:
                                        {
                                            prob_ = std::stod(word);
                                            break;
                                        }
                                    default:
                                        {
                                            break;
                                        }
                                }
                                ++id_counter;
                            }
                            if (running_marginal_probs[gene_][species_].find(count_key) !=
                                running_marginal_probs[gene_][species_].end()) {
                                running_marginal_probs[gene_][species_][count_key] +=
                                    prob_;
                            } else {
                                running_marginal_probs[gene_][species_][count_key] =
                                    prob_;
                            }
                        }
                        ++line_counter;

                    }

                    std::ofstream outfile;
                    outfile.open(filepath);
                    outfile << "count_key" << "," << "prob" << "\n";
                    for (auto &it : running_marginal_probs[gene_][species_]) {
                        // if ((it != running_marginal_probs[gene_][species_].end()) &&
                        //     (std::next(it) == running_marginal_probs[gene_][species_].end())) {
                        //     outfile << std::to_string(it->first) << "," <<
                        //         std::to_string(it->second);
                        // } else {
                        //     outfile << std::to_string(it->first) << "," <<
                        //         std::to_string(it->second) << "\n";
                        // }
                        outfile << std::to_string(it.first) << "," <<
                            std::to_string(it.second) << "\n";
                    }
                    running_marginal_probs[gene_][species_].clear();
                    outfile.close();

                } else {
                    std::ofstream outfile;
                    outfile.open(filepath);
                    outfile << "count_key" << "," << "prob" << "\n";

                    for (auto &it : running_marginal_probs[gene_][species_]) {
                        // if ((it != running_marginal_probs[gene_][species_].end()) &&
                        //     (std::next(it) == running_marginal_probs[gene_][species_].end())) {
                        //     outfile << std::to_string(it->first) << "," <<
                        //         std::to_string(it->second);
                        // } else {
                        //     outfile << std::to_string(it->first) << "," <<
                        //         std::to_string(it->second) << "\n";
                        // }
                        outfile << std::to_string(it.first) << "," <<
                            std::to_string(it.second) << "\n";
                    }
                    running_marginal_probs[gene_][species_].clear();
                    outfile.close();
                }


                for (unsigned int gene_2 = 0; gene_2 < num_genes; ++gene_2) {
                    for (int species_2 = 0; species_2 < reactions[gene_2].molecule_count_cur.size(); ++species_2) {
                        if ((gene_ == gene_2) && (species_ == species_2)) {

                        } else {
                            std::string filepath = filepath_probs + connector_ +
                                std::to_string(gene_) + connector_ +
                                std::to_string(species_) + connector_ +
                                    std::to_string(gene_2) + connector_ +
                                    std::to_string(species_2) + file_ext;

                            if (!if_first_write) {
                                std::ifstream prob_state_file(filepath);
                                std::string row_text;
                                std::string count_key;
                                double prob_;
                                std::string word;
                                unsigned int line_counter = 0;
                                while (std::getline(prob_state_file, row_text)) {
                                    if (line_counter > 0) {
                                        std::istringstream str_stream(row_text);

                                        unsigned int id_counter = 0;
                                        while (std::getline(str_stream, word, ',')) {
                                            switch(id_counter) {
                                                case 0:
                                                    {
                                                        count_key = word;
                                                        break;
                                                    }
                                                case 1:
                                                    {
                                                        prob_ = std::stod(word);
                                                        break;
                                                    }
                                                default:
                                                    {
                                                        break;
                                                    }
                                            }
                                            ++id_counter;
                                        }
                                        if (running_joint_probs[count_].find(count_key) !=
                                            running_joint_probs[count_].end()) {
                                            running_joint_probs[count_][count_key] +=
                                                prob_;
                                        } else {
                                            running_joint_probs[count_][count_key] = prob_;
                                        }
                                    }
                                    ++line_counter;

                                }

                                std::ofstream outfile;
                                outfile.open(filepath);
                                outfile << "count_key" << "," << "prob" << "\n";

                                for (auto &it : running_joint_probs[count_]) {
                                    // if ((it != running_joint_probs[count_].end()) &&
                                    //     (std::next(it) == running_joint_probs[count_].end())) {
                                    //     outfile << std::to_string(it->first) << "," <<
                                    //         std::to_string(it->second);
                                    // } else {
                                    //     outfile << std::to_string(it->first) << "," <<
                                    //         std::to_string(it->second) << "\n";
                                    // }
                                    outfile << it.first << "," <<
                                        std::to_string(it.second) << "\n";
                                }
                                running_joint_probs[count_].clear();
                                outfile.close();

                            } else {
                                std::ofstream outfile;
                                outfile.open(filepath);
                                outfile << "count_key" << "," << "prob" << "\n";

                                for (auto &it : running_joint_probs[count_]) {
                                    // if ((it != running_joint_probs[count_].end()) &&
                                    //     (std::next(it) == running_joint_probs[count_].end())) {
                                    //     outfile << std::to_string(it->first) << "," <<
                                    //         std::to_string(it->second);
                                    // } else {
                                    //     outfile << std::to_string(it->first) << "," <<
                                    //         std::to_string(it->second) << "\n";
                                    // }
                                    outfile << it.first << "," <<
                                        std::to_string(it.second) << "\n";
                                }
                                running_joint_probs[count_].clear();
                                outfile.close();
                            }
                            count_ += 1;
                        }
                    }
                }
            }
        }
    }

    void GillespieSSA::upate_running_probs (double total_time_, double next_time_step,
        unsigned int &history_num, bool &if_first_write) {
        unsigned int  count_ = 0;
        for (unsigned int gene_ = 0; gene_ < num_genes; ++gene_) {
            for (int species_ = 0; species_ < reactions[gene_].molecule_count_cur.size(); ++species_) {
                unsigned int prev_count = reactions[gene_].molecule_count_cur[species_];

                if (running_marginal_probs[gene_][species_].find(prev_count) !=
                    running_marginal_probs[gene_][species_].end()) {
                    running_marginal_probs[gene_][species_][prev_count] +=
                        next_time_step / total_time_;
                } else {
                    running_marginal_probs[gene_][species_][prev_count] =
                        next_time_step / total_time_;
                }

                for (unsigned int gene_2 = 0; gene_2 < num_genes; ++gene_2) {
                    for (int species_2 = 0; species_2 < reactions[gene_2].molecule_count_cur.size(); ++species_2) {
                        if ((gene_ == gene_2) && (species_ == species_2)) {

                        } else {
                            unsigned int prev_count_2 =
                                reactions[gene_2].molecule_count_cur[species_2];
                            std::string connector_ = "_";
                            std::string joint_counts = std::to_string(prev_count) +
                                connector_ + std::to_string(prev_count_2);

                            if (running_joint_probs[count_].find(joint_counts) !=
                                running_joint_probs[count_].end()) {
                                running_joint_probs[count_][joint_counts] +=
                                    next_time_step / total_time_;
                            } else {
                                running_joint_probs[count_][joint_counts] =
                                    next_time_step / total_time_;
                            }
                            count_ += 1;
                        }

                    }
                }
            }
        }

        if (history_num >= running_probs_buffer_size + 1) {
            output_probs_to_file(if_first_write);
            history_num = 0;
            if (if_first_write) {
                if_first_write = false;
            }
        }

    }

    void GillespieSSA::simulate (bool compute_statistics, std::string statistics_file_full,
        bool verbose, bool cell_cycle_sim_frozen) {
        // start_molecule_count_history_file();
        double count_rxns_num = 0;
        double statistics_start_time = 0;
        bool is_statistics_start_time_set = false;
        unsigned int which_random_time_saved = 0;
        is_steady_state_reached = false;
        unsigned int cur_rxn_count = 0;
        unsigned int rxn_count_after_burn_in = 0;
        double time_since_last_division = 0;
        std::string prev_cell_cycle_state = "G1";

        if (compute_statistics) {
            start_statistics_file(statistics_file_full);
            set_size_statistics_containers();
            // set_size_statistics_history_containers();
            statistics_start_time = 0;
            is_statistics_start_time_set = false;
        }else{
            is_steady_state_reached = true;
        }
        // time_history.clear();
        // time_history.push_back(0);
        std::fill(time_history.begin(), time_history.end(), 0);
        std::fill(rxn_history.begin(), rxn_history.end(), "transcription");

        update_burst_size_init();
        init_rxn_order();
        std::string tmp = get_cur_cell_cycle_state();
        // std::cout << "cur phase = " << reactions[0].rxn_rates["gene on"] << std::endl;
        init_cell_cycle_state (generator[0], 0);
        std::string cur_cell_cycle_phase = get_cur_cell_cycle_state();
        // std::cout << "cur phase2 = " << reactions[0].rxn_rates["gene on"] << std::endl;
        compute_total_propensity();
        set_count_rxns_fired(count_rxns, compute_statistics, stop_rxn_count);
        if (save_at_random_times && is_steady_state_reached) {
            find_random_times_to_save (generator[0],
                burn_in, max_time);
        }

        // std::random_device rd;
        // std::vector<std::uint_least32_t> rd_seeds = {rd(), rd(), rd(), rd()};
        // std::vector<std::uint_least32_t> rd_seeds =
        //     {random_seeds[0], random_seeds[1], random_seeds[2], random_seeds[3]};
        // std::vector<std::uint_least32_t> rd_seeds(random_seeds.size());
        // for (std::size_t b = 0; b < random_seeds.size(); ++b) {
        //     rd_seeds[b] = random_seeds[b];
        // }
        // std::seed_seq sd(rd_seeds.begin(), rd_seeds.end());
        // thread_local static RNG generator{sd};]
        bool reached_burn_in_rxn = false;
        bool stop_sim = false;
        bool reached_rxn_count = false;
        bool simulation_ended = false;
        int num_history = 0;
        int num_save_loop = 0;
        double total_time = 0;
        double cur_time = 0;
        double next_time_step = 0;
        std::string next_rxn_name = "transcription";
        double points_collected_prev = 0;
        cur_cell_cycle_phase = get_cur_cell_cycle_state();
        update_molecule_count_history(num_history, num_save_loop,
            simulation_ended, next_time_step, next_rxn_name, 
            time_since_last_division);


        // set params for running marginal and joint probs
        unsigned int history_num = 0;
        bool if_first_write = true;

        while (!stop_sim) {
            // std::cout << "start, max time = " << max_time << " burn in = " <<
            //     burn_in << std::endl;
            // cur_rxn_count += 1;
            if (compute_statistics && reached_burn_in_rxn) {
                rxn_count_after_burn_in += 1;
                if (!is_steady_state_reached && save_at_random_times) {
                    max_time = total_time + max_time - burn_in;
                    burn_in = total_time;
                    find_random_times_to_save (generator[0],
                        burn_in, max_time);
                    is_steady_state_reached = true;
                }
            }
            // if (reactions[0].molecule_count_cur[0] == 0 && reactions[0].propensity_vals["mRNA decay"] > 0) {
            //     std::cout << "error found0time " <<
            //     reactions[0].propensity_vals["mRNA decay"] << " " <<
            //     total_propensity << " time = " << total_time << std::endl;
            // }
            next_time_step = sample_time_step(generator[0]);
            if (reached_burn_in_rxn && compute_statistics) {
                if (!is_statistics_start_time_set) {
                    statistics_start_time = total_time;
                    is_statistics_start_time_set = true;
                }
                upate_running_statistics (total_time - statistics_start_time,
                    next_time_step);

                // upate_running_probs (total_time_norm, next_time_step,
                //     history_num, if_first_write);
                // ++history_num;
                // upate_running_statistics_without_time_step(count_rxns_num);
                count_rxns_num += 1;
            }
            // if (reactions[0].molecule_count_cur[0] == 0 && reactions[0].propensity_vals["mRNA decay"] > 0) {
            //     std::cout << "error found0time1 " <<
            //     reactions[0].propensity_vals["mRNA decay"] << " " <<
            //     total_propensity << std::endl;
            // }
            std::vector<std::string> GRN_out_changed;
            cur_time = total_time;
            total_time += next_time_step;
            if (verbose) {
                if (int(total_time / 1000) - int(cur_time / 1000) == 1) {
                    std::cout << "time = " << total_time << " saved point = " << which_random_time_saved << std::endl;
                }
            }
            // std::cout << "reached here 1 = " << std::endl;
            if (total_time < max_time || (count_rxns && !reached_rxn_count) ||
                (compute_statistics && !reached_rxn_count)) {
                // if (reactions[0].molecule_count_cur[0] == 0 && reactions[0].propensity_vals["mRNA decay"] > 0) {
                //     std::cout << "error found0 " <<
                //     reactions[0].propensity_vals["mRNA decay"] << " " <<
                //     total_propensity << std::endl;
                // }
                prev_cell_cycle_state = cur_cell_cycle_phase;
                update_cell_cycle_state(total_time, cur_time, generator[0]);
                cur_cell_cycle_phase = get_cur_cell_cycle_state();
                if (!cell_cycle_sim_frozen) {
                    if ((prev_cell_cycle_state == "G1") && 
                        (cur_cell_cycle_phase == "G1")) {
                        time_since_last_division += next_time_step;
                        }
                    if ((prev_cell_cycle_state == "G1") && 
                        (cur_cell_cycle_phase == "G2")) {
                        time_since_last_division += next_time_step;
                        }
                    if ((prev_cell_cycle_state == "G2") && 
                        (cur_cell_cycle_phase == "G1")) {
                        time_since_last_division = 0;
                        }
                    if ((prev_cell_cycle_state == "G2") && 
                        (cur_cell_cycle_phase == "G2")) {
                        time_since_last_division += next_time_step;
                        }
                }
                // if (reactions[0].molecule_count_cur[0] == 0 && reactions[0].propensity_vals["mRNA decay"] > 0) {
                //     std::cout << "error found01 " <<
                //     reactions[0].propensity_vals["mRNA decay"] << " " <<
                //     total_propensity << std::endl;
                // }

                int next_rxn = sample_next_rxn(generator[0]);
                std::string connector_ = ":";
                next_rxn_name =
                    std::to_string(rxn_order[next_rxn].gene_id) +
                    connector_ + rxn_order[next_rxn].rxn_name;
                update_burst_size (generator[0], next_rxn);
                reached_burn_in_rxn =
                    update_rxn_count (next_rxn, stop_sim, reached_rxn_count,
                        compute_statistics, burn_in_rxn_count);
                // std::cout << "reached here 2 = " << std::endl;
                GRN_out_changed = update_fired_reaction(next_rxn);
                update_dependent_count_propensity(next_rxn, GRN_out_changed);
                // if (reactions[0].molecule_count_cur[0] == 0 && reactions[0].propensity_vals["mRNA decay"] > 0) {
                //     std::cout << "error found3 " <<
                //     reactions[0].propensity_vals["mRNA decay"] << " " <<
                //     total_propensity << std::endl;
                // }
                // time_history.push_back(next_time_step);
                if (save_at_time_interval || save_at_random_times) {
                    if (save_at_time_interval) {
                        save_molecule_count_at_interval(cur_time, total_time,
                            points_collected_prev, next_rxn_name);
                    }else{
                        if (save_at_random_times && is_steady_state_reached) {
                            save_molecule_count_at_random_times(cur_time, total_time,
                                which_random_time_saved, next_rxn, 
                                time_since_last_division);
                        }
                    }
                }else{
                    update_molecule_count_history(num_history, num_save_loop,
                        simulation_ended, next_time_step, next_rxn_name, 
                        time_since_last_division);
                }
                // if (reactions[0].molecule_count_cur[0] == 0 && reactions[0].propensity_vals["mRNA decay"] > 0) {
                //     std::cout << "error foundlast " <<
                //     reactions[0].propensity_vals["mRNA decay"] << " " <<
                //     total_propensity << std::endl;
                // }
                // int rxn_found = 0;
                // for (int i = 0; i < rxn_order.size(); ++i) {
                //     if (rxn_order[i].rxn_name == "mRNA decay" &&
                //         rxn_order[i].gene_id == 0) {
                //             rxn_found = i;
                //         }
                // }
                // if (reactions[0].molecule_count_cur[0] == 0 &&
                //     reactions[0].propensity_vals["mRNA decay"] != 0) {
                //         std::cout << "aha!! 1 " << reactions[0].propensity_vals["mRNA decay"] << std::endl;
                //     }
                //
                // if (reactions[0].molecule_count_cur[0] == 0 &&
                //     rxn_order[rxn_found].propensity_val != 0) {
                //         std::cout << "aha!! 2 " << rxn_order[rxn_found].propensity_val << std::endl;
                //     }
            }else{
                // if (compute_statistics) {
                //     write_statistics_to_file(statistics_file_full);
                // }
                simulation_ended = true;
                // if (!save_at_time_interval) {
                //     update_molecule_count_history(num_history, num_save_loop,
                //         simulation_ended);
                // }
                stop_sim = true;
            }
        }
        if (compute_statistics) {
            write_statistics_to_file(statistics_file_full);

            // if (history_num > 1) {
            //     output_probs_to_file (if_first_write);
            // }
        }
        if (!save_at_time_interval && !save_at_random_times) {
            update_molecule_count_history(num_history, num_save_loop,
                simulation_ended, next_time_step, "transcription", 
                time_since_last_division);
        }
    }
}
