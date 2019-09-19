#!/bin/bash

canu -correct -p CanuCelegans -d CanuCelegans genomeSize=100286401 -pacbio-raw ../../longReads/simCelegans.fasta --stopOnReadQuality=false --corOutCoverage=300 --useGrid=false
cd CanuCelegans
gzip -d CanuCelegans.correctedReads.fasta.gz
cd ..
