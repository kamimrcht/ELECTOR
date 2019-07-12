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
cd ../miniasm;
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

cd ../src/split;
make;
cp masterSplitter ../../bin;
cp Donatello ../../bin;
