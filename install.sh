#!/bin/bash

mkdir bin;

cd src/simulator;
make;
cp simulator ../../bin;

cd ../poa-graph;
make poa;
cp poa ../../bin;

cd ../utils;
make;
cp fq2fa ../../bin;

cd ../../minimap2;
make;
cp minimap2 ../bin/
cd ../miniasm;
cp miniasm ../bin/
make;

cd ../htslib;
autoheader;
autoconf;
./configure --disable-lzma;
make;

cd ../samtools;
autoheader;
autoconf -Wno-syntax;
./configure --disable-lzma;
make;
cp samtools ../bin;

cd ../src/split;
make;
cp masterSplitter ../../bin;
cp Donatello ../../bin;
