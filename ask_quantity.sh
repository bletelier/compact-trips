#!/bin/bash

cd build/
for input in coruna-25 coruna-50 coruna-75 coruna-100
do
    for exes in baseline_ask_quantity PTR_ask_quantity PTR-C_ask_quantity TTR_ask_quantity TTR-C_ask_quantity
    do
        mkdir -p ../cds/$input
        echo "Executing ask_quantiy of $exes using $input.txt"
        ./$exes $input
        echo ""
    done
done
cd ..
