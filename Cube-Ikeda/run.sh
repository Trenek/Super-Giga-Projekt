#!/usr/bin/env bash
set -e

# compile C++
g++ main.cpp -o cubicIkeda  `capd-config --cflags --libs`

# run C++ → generate CSV
./cubicIkeda output/result.csv

# run Python → generate plot
python3 plot_traj.py