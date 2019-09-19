#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected ../../longReads/simEcoli -corrected ../../correction/Daccord/DaccordEcoli.fasta -reference ../../references/Ecoli.fasta -simulator simlord -output DaccordEcoli -split -dazzDb ../../correction/Daccord/EcoliDb.db -corrector daccord
cd reproduce_manuscript_results/evaluation/Daccord/
rm corrected_* reference_* uncorrected_*
