#!/bin/bash
make clean
make cm
./cm count -C -k 22 -h 2 -w 4 -r 5 -t 3 -o output.bin -fa rymv.sim.fa
./cm query -f output.bin -q query.txt -o query_result.txt