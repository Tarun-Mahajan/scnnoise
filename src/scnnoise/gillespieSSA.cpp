#include <vector>
#include <algorithm>
#include <iostream>
#include <omp.h>
#include <string>
#include <random>
#include <math.h>
#include "graph.hpp"
#include "scnnoise.hpp"
#include "GillespieSSA.hpp"

namespace ScnnoiseInterface {
  namespace SimulatorGillespieSSA {
    /* function definitions */
    // Constructor
    GillespieSSA::GillespieSSA (int num_rxns, int num_genes,
      const std::vector<int> num_species_gene_type,
      const std::vector<int> num_rxns_gene_type, double max_time,
      bool save_timeseries, int num_timepoints_save,
      std::string count_save_file):
      scnnoise (num_rxns, num_genes, num_species_gene_type,
        num_rxns_gene_type, max_time, save_timeseries, num_timepoints_save,
        count_save_file) {
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

      for (int i = 0; i < num_rxns; ++i) {
        selector -= rxn_order[i].propensity_val;
        if (selector <= 0) {
          rxn_selected = i;
          break;
        }
      }
      sort_reaction(rxn_selected);
      // if (rxn_selected > 0) {
      //   std::swap(rxn_order[rxn_selected - 1], rxn_order[rxn_selected]);
      //   rxn_selected -= 1;
      // }
      return rxn_selected;
    }

    void GillespieSSA::update_fired_reaction_reactants(int gene_selected,
      int rxn_index, int &count_not_changed_reactants,
      std::vector<bool> &flag_changed_product_count,
      std::vector<bool> &GRN_out_changed,
      std::vector<int> &rxn_selected_reactants,
      std::vector<int> &rxn_selected_products,
      std::vector<int> &rxn_selected_reactants_stoichio,
      std::vector<int> &rxn_selected_products_stoichio) {
        for (auto r = rxn_selected_reactants.begin(); r != rxn_selected_reactants.end(); ++r) {
          int reactant_index = std::distance(rxn_selected_reactants.begin(), r);
          std::vector<int>::iterator it =
          std::find(rxn_selected_products.begin(), rxn_selected_products.end(),
          *r);
          if (it != rxn_selected_products.end()) {
            int product_index = std::distance(rxn_selected_products.begin(), it);
            flag_changed_product_count[product_index] = true;
            int count_change = rxn_selected_products_stoichio[product_index] -
              rxn_selected_reactants_stoichio[reactant_index]
            if (count_change) {
              reactions[gene_selected].molecule_count_cur[*r] += count_change;
              std::vector<int>::iterator it1 =
              std::find(reactions[gene_selected].GRN_rxn_OUT.begin(),
              reactions[gene_selected].GRN_rxn_OUT.end(),
              *r);
              if (it1 != reactions[gene_selected].GRN_rxn_OUT.end()) {
                int out_index = std::distance(reactions[gene_selected].GRN_rxn_OUT.begin(), it1);
                GRN_out_changed[out_index] = true;
              }
              // if (*r == reactions[gene_selected].GRN_rxn_OUT) {
              //   GRN_out_changed = true;
              // }
            }else{
              count_not_changed_reactants += 1;
            }
          }else{
            reactions[gene_selected].molecule_count_cur[*r] -=
              rxn_selected_reactants_stoichio[reactant_index];
            std::vector<int>::iterator it1 =
            std::find(reactions[gene_selected].GRN_rxn_OUT.begin(),
            reactions[gene_selected].GRN_rxn_OUT.end(),
            *r);
            if (it1 != reactions[gene_selected].GRN_rxn_OUT.end()) {
              int out_index = std::distance(reactions[gene_selected].GRN_rxn_OUT.begin(), it1);
              GRN_out_changed[out_index] = true;
            }
          }
        }
      }

    std::vector<bool> GillespieSSA::update_fired_Reaction (int rxn_selected) {
      /*
      Find the gene and reaction type for the selected reaction channel. Extract the reactants,
      products and their stoichiometric coefficients for the selected reaction channel.
      */
      int gene_selected = rxn_order[rxn_selected].gene_id;
      // int rxn_type_selected = rxn_order[rxn_selected].rxn_type;
      int rxn_index = rxn_order[rxn_selected].rxn_type;
      // std::vector<int> rxn_temp = reactions[gene_selected].rxn_type;
      // std::vector<int>::iterator it = std::find(rxn_temp.begin(), rxn_temp.end(),
      //   rxn_type_selected);
      // int rxn_index = std::distance(rxn_temp.begin(), it);
      std::vector<int> rxn_selected_reactants =
        reactions[gene_selected].rxns[rxn_index].reactants;
      std::vector<int> rxn_selected_products =
        reactions[gene_selected].rxns[rxn_index].products;
      std::vector<int> rxn_selected_reactants_stoichio =
        reactions[gene_selected].rxns[rxn_index].reactants_stoichio;
      std::vector<int> rxn_selected_products_stoichio =
        reactions[gene_selected].rxns[rxn_index].products_stoichio;


      std::vector<bool> GRN_out_changed(reactions[gene_selected].GRN_rxn_OUT.size(), false);
      /*
      Subtract current reaction propensity from total reaction propensity for
      computing updated total propensity later.
      */
      total_propensity -= reactions[gene_selected].rxns[rxn_index].propensity_val;

      /*
      Update molecular count for reactants which change by firing the selected
      reaction channel.
      */
      int count_not_changed_reactants = 0;
      std::vector<bool> flag_changed_product_count(rxn_selected_products.size(), false);
      update_fired_reaction_reactants(gene_selected, rxn_index,
        count_not_changed_reactants, flag_changed_product_count,
        GRN_out_changed, rxn_selected_reactants, rxn_selected_products,
        rxn_selected_reactants_stoichio, rxn_selected_products_stoichio);

      /*
      Update molecular count for products which change by firing the selected
      reaction channel.
      */
      for (auto r = rxn_selected_products.begin(); r != rxn_selected_products.end(); ++r) {
        int product_index = std::distance(rxn_selected_products.begin(), r);
        if (!flag_changed_product_count[product_index]) {
          reactions[gene_selected].molecule_count_cur[*r] +=
            rxn_selected_products_stoichio[product_index];
            // if (*r == reactions[gene_selected].GRN_rxn_OUT) {
            //   GRN_out_changed = true;
            // }

            std::vector<int>::iterator it1 =
            std::find(reactions[gene_selected].GRN_rxn_OUT.begin(),
            reactions[gene_selected].GRN_rxn_OUT.end(),
            *r);
            if (it1 != reactions[gene_selected].GRN_rxn_OUT.end()) {
              int out_index = std::distance(reactions[gene_selected].GRN_rxn_OUT.begin(), it1);
              GRN_out_changed[out_index] = true;
            }
        }
      }

      /*
      Update propensity for the selected reaction channel.
      */
      if (count_not_changed_reactants < rxn_selected_reactants.size()) {
        double new_propensity = reactions[gene_selected].rxns[rxn_index].rxn_rate;
        for (auto r = rxn_selected_reactants.begin(); r != rxn_selected_reactants.end(); ++r) {
          int reactant_index = std::distance(rxn_selected_reactants.begin(), r);
          int reactant_stoichio = rxn_selected_reactants_stoichio[reactant_index];
          new_propensity *= (factorial(reactions[gene_selected].molecule_count_cur[*r])/
            factorial(reactions[gene_selected].molecule_count_cur[*r] - reactant_stoichio));
        }
        std::vector<int>::iterator it1 =
        std::find(reactions[gene_selected].GRN_rxn_IN.begin(),
        reactions[gene_selected].GRN_rxn_IN.end(),
        rxn_index);
        if (it1 != reactions[gene_selected].GRN_rxn_IN.end()) {
          new_propensity *= regulation_function(gene_selected, rxn_index);
        }


        // if (rxn_index == reactions[gene_selected].GRN_rxn_IN) {
        //   new_propensity *= regulation_function(gene_selected);
        // }
        reactions[gene_selected].rxns[rxn_index].propensity_val = new_propensity;
        rxn_order[rxn_selected].propensity_val = new_propensity;

      }
      total_propensity += reactions[gene_selected].rxns[rxn_index].propensity_val;
      return GRN_out_changed;
    }

    void GillespieSSA::update_dependent_count_propensity (int rxn_selected,
      std::vector<bool> &GRN_out_changed) {

        /*
        Update propensity for dependent reactions belonging to the same gene as the
        reaction selected for firing
        */
        int gene_selected = rxn_order[rxn_selected].gene_id;
        int rxn_index = rxn_order[rxn_selected].rxn_type;
        int gene_type = reactions[gene_selected].gene_type;
        vector<int> rxn_selected_children;
        rxn_selected_children.reserve(gene_rxn_dependency[gene_type].get_size());
        gene_rxn_dependency[gene_type].find_children(rxn_index, rxn_selected_children);
        std::vector<int> rxn_selected_reactants;
        rxn_selected_reactants.reserve(reactions[gene_selected].molecule_count_cur.size());
        std::vector<int> rxn_selected_reactants_stoichio;
        rxn_selected_reactants_stoichio.reserve(reactions[gene_selected].molecule_count_cur.size());

        for (auto &rxn : rxn_selected_children) {
          total_propensity -= reactions[gene_selected].rxns[rxn].propensity_val;
          rxn_selected_reactants =
            reactions[gene_selected].rxns[rxn].reactants;
          rxn_selected_reactants_stoichio =
            reactions[gene_selected].rxns[rxn].reactants_stoichio;

          double new_propensity = reactions[gene_selected].rxns[rxn].rxn_rate;
          for (auto r = rxn_selected_reactants.begin(); r != rxn_selected_reactants.end(); ++r) {
            int reactant_index = std::distance(rxn_selected_reactants.begin(), r);
            int reactant_stoichio = rxn_selected_reactants_stoichio[reactant_index];
            new_propensity *= (factorial(reactions[gene_selected].molecule_count_cur[*r])/
              factorial(reactions[gene_selected].molecule_count_cur[*r] - reactant_stoichio));
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
              if (out_species == reactions[gene_selected].GRN_rxn_OUT[g_index]) {
                total_propensity -= reactions[gene_children[g]].rxns[rxn].propensity_val;
                rxn_selected_reactants =
                  reactions[gene_selected[g]].rxns[rxn].reactants;
                rxn_selected_reactants_stoichio =
                  reactions[gene_selected[g]].rxns[rxn].reactants_stoichio;

                double new_propensity = reactions[gene_children[g]].rxns[rxn].rxn_rate;
                for (auto r = rxn_selected_reactants.begin(); r != rxn_selected_reactants.end(); ++r) {
                  int reactant_index = std::distance(rxn_selected_reactants.begin(), r);
                  int reactant_stoichio = rxn_selected_reactants_stoichio[reactant_index];
                  new_propensity *= (factorial(reactions[gene_selected[g]].molecule_count_cur[*r])/
                    factorial(reactions[gene_selected[g]].molecule_count_cur[*r] - reactant_stoichio));
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
              }

            }
          }
        }
    }

    void GillespieSSA::update_molecule_count_history (int &num_history, int &num_save_loop) {
      if (num_history < num_timepoints_save) {
      }else{
        num_history = 0;
        num_save_loop += 1;
        if (save_timeseries) {
          std::ofstream outfile;
          outfile.open(count_save_file, std::ios_base::app);
          for (int id_time = 0; id_time < num_timepoints_save; ++id_time) {
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
      }
      for (int gene; gene < num_genes; ++gene) {
        for (int species; species < num_species_gene_type[gene]; ++species) {
          molecule_count_history[gene][species][num_history] =
            reactions[gene].molecule_count_cur[species];
        }
      }
    }

    double GillespieSSA::simulate () override {
      compute_total_propensity();

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
        bool GRN_out_changed = false;
        double total_time = std::accumulate(time_history.begin(), time_history.end(),
          decltype(time_history)::value_type(0)) + next_time_step;
        if (total_time < max_time) {
          update_fired_Reaction(next_rxn, GRN_out_changed);
          update_dependent_count_propensity(next_rxn, GRN_out_changed);
          time_history.push_back(next_time_step);
          update_molecule_count_history(num_history, num_save_loop);
        }else{
          stop_sim = true;
        }
      }
    }
  }
}
