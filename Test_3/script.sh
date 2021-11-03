#!/bin/bash
rm cm
g++ main.cpp -o cm
./cm count -C -k 22 -h 1000 -w 10000 -r 5 -t 3 -o output.bin -fa rymv.sim.fa
./cm query -f output.bin -q query.txt -o query_result.txt