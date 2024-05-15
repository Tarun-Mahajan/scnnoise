#include <iostream>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <sstream>
// #include <string>
// #include "../src/scnnoise/graph_derived.hpp"
#include "utils.hpp"

void create_GRN_from_file (GraphSpace::GRN &gene_net, std::string filepath,
    std::map<std::string, int> gene_rev_map) {
    std::ifstream GRN_file(filepath);
    std::string row_text;
    std::vector<int> GRN_int_param;
    std::vector<double> GRN_params;
    std::vector<bool> GRN_activation;
    std::string word;
    while (std::getline(GRN_file, row_text)) {
        GRN_int_param.clear();
        GRN_params.clear();
        GRN_activation.clear();
        std::istringstream str_stream(row_text);

        unsigned int id_counter = 0;
        while (std::getline(str_stream, word, ',')) {
            if (id_counter < 2 || (id_counter >=5 && id_counter <= 6)) {
                if (id_counter < 2) {
                    GRN_int_param.push_back(gene_rev_map[word]);
                }else{
                    GRN_int_param.push_back(std::stoi(word));
                }
            }else{
                if (id_counter == 7) {
                    std::transform(word.begin(), word.end(), word.begin(), ::tolower);
                    std::istringstream is(word);
                    bool bool_val;
                    is >> std::boolalpha >> bool_val;
                    GRN_activation.push_back(bool_val);
                }else{
                    GRN_params.push_back(std::stod(word));
                }
            }
            id_counter += 1;
        }
        gene_net.add_edge_kinetics(GRN_int_param[0], GRN_int_param[1],
            GRN_params[0], GRN_params[1], GRN_params[2],
            GRN_int_param[2], GRN_int_param[3], GRN_activation[0]);

    }
}
