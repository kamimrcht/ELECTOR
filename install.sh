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
mv fq2fa ../../bin;
cd ../../minimap2;
make;
cd ../miniasm;
make;
cd ../bwa;
make;
cd ../htslib;
autoheader;
autoconf;
./configure;
make;
cd ../samtools;
autoheader:
autoconf -Wno-syntax;
./configure;
make;
cd ../src/split;
make;
mv masterSplitter ../../bin;
mv Donatello ../../bin;
