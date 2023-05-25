#!/bin/bash

cd build/

for exes in k3-v2_create
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
