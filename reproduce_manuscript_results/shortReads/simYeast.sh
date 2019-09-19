#!/bin/bash

art_illumina -ss MSv3 -i ../references/Yeast.fasta -l 250 -f 50 -o simYeast
./../../bin/fq2fa simYeast.fq > simYeast.fa
