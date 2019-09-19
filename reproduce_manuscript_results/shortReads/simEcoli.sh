#!/bin/bash

art_illumina -ss MSv3 -i ../references/Ecoli.fasta -l 250 -f 50 -o simEcoli
./../../bin/fq2fa simEcoli.fq > simEcoli.fa
