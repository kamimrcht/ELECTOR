#!/bin/bash

./ConvertToPacBio_q2a.py ../../longReads/simYeast.fasta
mv LR.fasta simYeast.fasta
fasta2DB YeastDb simYeast.fasta
daligner YeastDb YeastDb
daccord -t$(nproc) YeastDb.YeastDb.las YeastDb.db > DaccordYeast.fasta
rm simYeast.fasta
