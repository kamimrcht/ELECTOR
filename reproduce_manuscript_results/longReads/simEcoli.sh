#!/bin/bash

source activate simlord
simlord --read-reference ../references/Ecoli.fasta -t 0.305 -pi 0.22 -pd 0.08 -ps 0.01 --coverage 20 --sam-output simEcoli.sam simEcoli
./../../bin/fq2fa simEcoli.fastq > simEcoli.fasta
source deactivate
