#!/bin/bash

cd build/
for input in coruna-25 coruna-50 coruna-75 coruna-100
do
    for exes in baseline_ask_quantity nok3-v1_ask_quantity nok3-v2_ask_quantity k3-v1_ask_quantity k3-v2_ask_quantity
    do
        mkdir -p ../cds/$input
        echo "Executing ask_quantiy of $exes using $input.txt"
        ./$exes $input
        echo ""
    done
done
cd ..
