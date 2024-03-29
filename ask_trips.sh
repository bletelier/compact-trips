#!/bin/bash

cd build/
for input in coruna-25 coruna-50 coruna-75 coruna-100
do
    for exes in baseline_ask_trips PTR_ask_trips PTR-C_ask_trips TTR_ask_trips TTR-C_ask_trips
    do
        mkdir -p ../cds/$input
        echo "Executing ask_trips of $exes using $input.txt"
        ./$exes $input
        echo ""
    done
done
cd ..
