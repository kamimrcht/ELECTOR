#!/bin/bash

./minia -in ../../shortReads/simYeast.fa
python2.7 runHALC.py -o ../../shortReads/simYeast.fa -t $(nproc) ../../longReads/simYeast.fasta simYeast.contigs.fa
rm simYeast.*
mv output HALCYeast
rm -rf temp
