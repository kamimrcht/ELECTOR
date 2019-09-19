#!/bin/bash

lordec-correct -T $(nproc) -i ../../longReads/simCelegans.fasta -2 ../../shortReads/Celegans.fa -k 23 -o LordecCelegans.fasta -s 3
lordec-trim-split -i LordecCelegans.fasta -o LordecCelegans.split.fasta
