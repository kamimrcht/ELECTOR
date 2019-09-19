#!/bin/bash

./minia -in ../../shortReads/Celegans.fa
python2.7 runHALC.py -o ../../shortReads/Celegans.fa -t $(nproc) ../../longReads/simCelegans.fasta Celegans.contigs.fa
rm Celegans.*
mv output HALCCelegans
rm -rf temp
