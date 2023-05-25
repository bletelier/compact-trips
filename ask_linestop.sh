#!/bin/bash

cd build/

for exes in baseline_ask_linestop PTR_ask_linestop PTR-C_ask_linestop TTR_ask_linestop TTR-C_ask_linestop
do
    for input in coruna-25 coruna-50 coruna-75 coruna-100
    do
        mkdir -p ../cds/$input
        echo "Executing ask_linestop of $exes using $input.txt"
        ./$exes $input
        echo ""
    done
done
cd ..
