#!/bin/bash

./minia -in ../../shortReads/Yeast.fa
python2.7 runHALC.py -o ../../shortReads/Yeast.fa -t $(nproc) ../../longReads/simYeast.fasta Yeast.contigs.fa
rm Yeast.*
mv output HALCYeast
rm -rf temp
