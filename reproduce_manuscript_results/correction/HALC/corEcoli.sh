#!/bin/bash

./minia -in ../../shortReads/simEcoli.fa
python2.7 runHALC.py -o ../../shortReads/simEcoli.fa -t $(nproc) ../../longReads/simEcoli.fasta simEcoli.contigs.fa
rm simEcoli.*
mv output HALCEcoli
rm -rf temp
