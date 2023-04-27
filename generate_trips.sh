#!/bin/bash

cd build/

for input in coruna-100 coruna-75 coruna-50 coruna-25
do
    echo "Executing $exes using trips/$input.txt"
    ./trip_generator $input < ../maps/$input.txt
    echo ""
done
cd ..
