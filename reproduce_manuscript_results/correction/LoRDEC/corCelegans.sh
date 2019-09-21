#!/bin/bash

lordec-correct -T $(nproc) -i ../../longReads/simCelegans.fasta -2 ../../shortReads/simCelegans.fa -k 23 -o LordecCelegans.fasta -s 3
lordec-trim-split -i LordecCelegans.fasta -o LordecCelegans.split.fasta
./../removeLineBreaks.sh LordecCelegans.split.fasta
