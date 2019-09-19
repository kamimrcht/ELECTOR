#!/bin/bash

mecat2pw -j 0 -d ../../longReads/simCelegans.fasta -o candidates.txt -w . -t $(nproc) -x 1
mecat2cns -i 0 -t $(nproc) -x 1 candidates.txt ../../longReads/simCelegans.fasta MECATCelegans.fasta
rm candidates.txt* fileindex.txt r_0 vol0
