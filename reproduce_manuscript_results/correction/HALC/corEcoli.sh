#!/bin/bash

./minia -in ../../shortReads/Ecoli.fa
python2.7 runHALC.py -o ../../shortReads/Ecoli.fa -t $(nproc) ../../longReads/simEcoli.fasta Ecoli.contigs.fa
rm Ecoli.*
mv output HALCEcoli
rm -rf temp
