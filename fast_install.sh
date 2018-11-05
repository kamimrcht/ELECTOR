#!/bin/bash

mkdir bin;

cd src/simulator;
make -j;
cp simulator ../../bin;

cd ../poa-graph;
make poa -j;
cp poa ../../bin;

cd ../utils;
make -j;
cp fq2fa ../../bin;

cd ../../minimap2;
make -j;
cd ../miniasm;
make -j;

cd ../bwa;
make -j;

cd ../htslib;
autoheader;
autoconf;
./configure --disable-lzma;
make -j;

cd ../samtools;
autoheader;
autoconf -Wno-syntax;
./configure --disable-lzma;
make -j;

cd ../src/split;
make -j;
cp masterSplitter ../../bin;
cp Donatello ../../bin;
