#!/bin/bash

lordec-correct -T $(nproc) -i ../../longReads/simYeast.fasta -2 ../../shortReads/simYeast.fa -k 19 -o LordecYeast.fasta -s 3
lordec-trim-split -i LordecYeast.fasta -o LordecYeast.split.fasta
