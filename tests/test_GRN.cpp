#include <catch2/catch_all.hpp>
#include <iostream>
#include <vector>
#include "test_utils.hpp"
#include "test_GRN_global_variables.hpp"
// #include "../src/scnnoise/graph_derived.hpp"

extern std::string GRN_filepath;
extern unsigned int num_nodes;

TEST_CASE("Test GRN class")
{
    std::cout << "Testing GRN class for network with " << num_nodes <<
        " genes" << std::endl;

    GraphSpace::GRN gene_net(num_nodes);
    create_GRN_from_file (gene_net, GRN_filepath);

    SECTION ("Test if node count is correct") {
        bool is_node_count_set = test_node_count (gene_net, num_nodes);
        REQUIRE(is_node_count_set == true);
    }

    SECTION ("Test if edge count is correct") {
        bool is_edge_count_set = test_edges_count (gene_net, GRN_filepath);
        REQUIRE(is_edge_count_set == true);
    }

    SECTION ("Test if edge properties are set") {
        std::vector<bool> are_edge_properties_set =
            test_edge_properties (gene_net, GRN_filepath);
        std::vector<bool> edge_properties_true(3, true);
        REQUIRE(are_edge_properties_set == edge_properties_true);
    }

    SECTION ("Test if function 'find_children' works") {
        bool is_find_children_working = test_find_children (gene_net);
        REQUIRE(is_find_children_working == true);
    }
}
