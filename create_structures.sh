#!/bin/bash

cd build/

for exes in baseline_create PTR_create PTR-C_create TTR_create TTR-C_create
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
