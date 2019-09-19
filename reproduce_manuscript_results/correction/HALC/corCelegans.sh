#!/bin/bash

./minia -in ../../shortReads/simCelegans.fa
python2.7 runHALC.py -o ../../shortReads/simCelegans.fa -t $(nproc) ../../longReads/simCelegans.fasta simCelegans.contigs.fa
rm simCelegans.*
mv output HALCCelegans
rm -rf temp
