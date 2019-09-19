#!/bin/bash

canu -correct -p CanuEcoli -d CanuEcoli genomeSize=4641652 -pacbio-raw ../../longReads/simEcoli.fasta --stopOnReadQuality=false --corOutCoverage=300 --useGrid=false
cd CanuEcoli
gzip -d CanuEcoli.correctedReads.fasta.gz
cd ..
