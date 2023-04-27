#!/bin/bash

cd build/
for input in coruna-100 coruna-75 coruna-50 coruna-25
do
    for exes in baseline_ask_quantity k3-v1_ask_quantity k3-v2_ask_quantity nok3-v1_ask_quantity nok3-v2_ask_quantity nok3-v3_ask_quantity
    do
        mkdir -p ../cds/$input
        echo "Executing ask_quantiy of $exes using $input.txt"
        ./$exes $input
        echo ""
    done
done
cd ..
