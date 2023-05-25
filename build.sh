#!/bin/bash
mkdir -p build
cd src/

for folder in baseline TTR-C TTR PTR-C PTR sintetic_trips_generator
do
    cd $folder
    make
    cd ..
done
