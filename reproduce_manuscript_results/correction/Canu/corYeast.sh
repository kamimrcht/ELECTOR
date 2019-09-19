#!/bin/bash

canu -correct -p CanuYeast -d CanuYeast genomeSize=12370681 -pacbio-raw ../../longReads/simYeast.fasta --stopOnReadQuality=false --corOutCoverage=300 --useGrid=false
cd CanuYeast
gzip -d CanuYeast.correctedReads.fasta.gz
cd ..
