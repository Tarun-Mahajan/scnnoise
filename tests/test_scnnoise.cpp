#include <catch.hpp>
#include "scnnoise.hpp"
#include "gillespieSSA.hpp"
#include "gillespieSDM.hpp"

TEST_CASE("Test scNNoiSE", )
{
    int num_rxns = 14;
    int num_genes = 2;
    const std::vector<int> num_species_gene_type(1, 5);
    const std::vector<int> num_rxns_gene_type(1, 7);
    double max_time = 10000;
    bool save_timeseries = true;
    int num_timepoints_save = 100;
    std::string count_save_file = "data/testing_temp.two_gene_test.csv";
    ScnnoiseInterface::GillespieSDM obj(num_rxns, num_genes,
                                       num_species_gene_type, num_rxns_gene_type,
                                       max_time, save_timeseries,
                                       num_timepoints_save, count_save_file);
    REQUIRE(obj.num_rxns == num_rxns);
}
