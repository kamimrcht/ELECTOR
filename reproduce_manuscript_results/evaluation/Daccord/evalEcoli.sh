#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected reproduce_manuscript_results/longReads/simEcoli -corrected reproduce_manuscript_results/correction/Daccord/DaccordEcoli.fasta -reference reproduce_manuscript_results/references/Ecoli.fasta -simulator simlord -output DaccordEcoli -split -dazzDb reproduce_manuscript_results/correction/Daccord/EcoliDb.db -corrector daccord
cd reproduce_manuscript_results/evaluation/Daccord/
rm corrected_* reference_* uncorrected_*
