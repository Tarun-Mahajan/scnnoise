#include <catch2/catch_all.hpp>
#include "scnnoise.hpp"
#include "gillespieSSA.hpp"
#include "gillespieSDMnoCellCycle.hpp"

TEST_CASE("Test scNNoiSE")
{
    int num_rxns = 14;
    int num_genes = 2;
    std::vector<int> num_species_gene_type(1, 5);
    std::vector<int> num_rxns_gene_type(1, 7);
    double max_time = 10000;
    bool save_timeseries = true;
    int num_timepoints_save = 100;
    std::string count_save_file = "data/testing_temp.two_gene_test.csv";
    ScnnoiseInterface::gillespieSDMnoCellCycle obj(num_rxns, num_genes,
                                       num_species_gene_type, num_rxns_gene_type,
                                       max_time, save_timeseries,
                                       num_timepoints_save, count_save_file);
    REQUIRE(num_rxns == num_rxns);
    REQUIRE(obj.num_genes == num_genes);
    REQUIRE(obj.num_species_gene_type == num_species_gene_type);
    REQUIRE(obj.num_rxns_gene_type == num_rxns_gene_type);
    REQUIRE(obj.max_time == max_time);
    REQUIRE(obj.save_timeseries == save_timeseries);
    REQUIRE(obj.num_timepoints_save == num_timepoints_save);
    REQUIRE(obj.count_save_file == count_save_file);
}
