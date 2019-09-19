#!/bin/bash

./ConvertToPacBio_q2a.py ../../longReads/simEcoli.fasta
mv LR.fasta simEcoli.fasta
fasta2DB EcoliDb simEcoli.fasta
daligner EcoliDb EcoliDb
daccord -t$(nproc) EcoliDb.EcoliDb.las EcoliDb.db > DaccordEcoli.fasta
rm simEcoli.fasta
./../removeLineBreaks.sh DaccordEcoli.fasta
