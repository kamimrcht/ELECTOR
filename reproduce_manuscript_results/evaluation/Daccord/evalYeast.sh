#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected ../../longReads/simYeast -corrected ../../correction/Daccord/DaccordYeast.fasta -reference ../../references/Yeast.fasta -simulator simlord -output DaccordYeast -split -dazzDb ../../correction/Daccord/YeastDb.db -corrector daccord
cd reproduce_manuscript_results/evaluation/Daccord/
rm corrected_* reference_* uncorrected_*
