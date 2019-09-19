#!/bin/bash

source activate simlord
simlord --read-reference ../references/Celegans.fasta -t 0.305 -pi 0.22 -pd 0.08 -ps 0.01 --coverage 20 --sam-output simCelegans.fastq.sam simCelegans
./../../bin/fq2fa simCelegans.fastq > simCelegans.fasta
source deactivate
