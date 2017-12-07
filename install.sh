#!/bin/bash

mkdir bin;
cd src/simulator;
make;
mv simulator ../../bin;
cd ../poa-graph;
make poa;
mv poa ../../bin;
cd ../utils;
make;
mv fa2fq ../../bin;
