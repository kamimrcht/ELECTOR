#!/bin/bash

mecat2pw -j 0 -d ../../longReads/simEcoli.fasta -o candidates.txt -w . -t $(nproc) -x 1
mecat2cns -i 0 -t $(nproc) -x 1 candidates.txt ../../longReads/simEcoli.fasta MECATEcoli.fasta
rm candidates.txt* fileindex.txt r_0 vol0
