#!/usr/bin/env bash
cd ../build
cmake ../tests
cmake --build . --target scnnoise_graph_test

GRN_filepath="../data/test_GRN/GRN_size_6_net_1.csv"
num_genes=6
./scnnoise_graph_test -g $GRN_filepath --num_nodes $num_genes
