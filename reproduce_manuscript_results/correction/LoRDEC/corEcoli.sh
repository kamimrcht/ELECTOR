#!/bin/bash

lordec-correct -T $(nproc) -i ../../longReads/simEcoli.fasta -2 ../../shortReads/simEcoli.fa -k 19 -o LordecEcoli.fasta -s 3
lordec-trim-split -i LordecEcoli.fasta -o LordecEcoli.split.fasta
./../removeLineBreaks.sh LordecEcoli.split.fasta
