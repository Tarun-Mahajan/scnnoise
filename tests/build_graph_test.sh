#!/usr/bin/env bash
cd ../build
cmake ../tests
cmake --build . --target scnnoise_graph_test
