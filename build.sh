#!/bin/bash
mkdir -p build
cd src/

for folder in baseline k3_v1 k3_v2 nok3_v1 nok3_v2 nok3_v3 sintetic_trips_generator
do
    cd $folder
    make
    cd ..
done
