#!/bin/bash

source activate simlord
simlord --read-reference ../references/Yeast.fasta -t 0.305 -pi 0.22 -pd 0.08 -ps 0.01 --coverage 20 --sam-output simYeast.fastq.sam simYeast
./../../bin/fq2fa simYeast.fastq > simYeast.fasta
source deactivate
