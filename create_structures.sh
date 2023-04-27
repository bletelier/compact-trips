#!/bin/bash

cd build/

for exes in baseline_create k3-v1_create k3-v2_create nok3-v1_create nok3-v2_create nok3-v3_create
do
    for input in coruna-100 coruna-75 coruna-50 coruna-25
    do
        mkdir -p ../cds/$input
        echo "Executing $exes using trips/$input.txt"
        ./$exes $input < ../trips/$input.txt
        echo ""
    done
done
cd ..
