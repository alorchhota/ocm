#!/bin/bash
rm ocm
g++ -std=c++17 -o ocm main.cpp
./ocm count -k 22 -h 7 -w 1048576 -n 4 -t 3 -o output.bin -fa rymv.sim.fa
./ocm query -f output.bin -q query.txt -o query_result.txt
