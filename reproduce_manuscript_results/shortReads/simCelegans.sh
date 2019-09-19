#!/bin/bash

art_illumina -ss MSv3 -i ../references/Celegans.fasta -l 250 -f 50 -o simCelegans
./../../bin/fq2fa simCelegans.fq > simCelegans.fa
