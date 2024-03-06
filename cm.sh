#!/bin/bash
rm cm
g++ -std=c++17 -o cm cm.cpp
#./cm count [-c = conservative update] [-r = do not canonicalize] -k <kmer-length> [-h <counters-height>] [-w <counters-width>] -o <out-ocm-file> -fa <input-fasta-file>
./cm count -c -k 22 -h 7 -w 1048576 -o output.bin -fa input/lhg22L20MC5x.fa
./cm query -f output.bin -q input/test_exact_count_lhg22L20MC5x_20000.txt -o output/query_result.csv
