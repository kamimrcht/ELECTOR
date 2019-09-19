#!/bin/bash

./ConvertToPacBio_q2a.py ../../longReads/simCelegans.fasta
mv LR.fasta simCelegans.fasta
fasta2DB CelegansDb simCelegans.fasta
daligner CelegansDb CelegansDb
daccord -t$(nproc) CelegansDb.CelegansDb.las CelegansDb.db > DaccordCelegans.fasta
rm simCelegans.fasta
./../removeLineBreaks.sh DaccordCelegans.fasta
