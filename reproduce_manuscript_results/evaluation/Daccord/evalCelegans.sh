#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected reproduce_manuscript_results/longReads/simCelegans -corrected reproduce_manuscript_results/correction/Daccord/DaccordCelegans.fasta -reference reproduce_manuscript_results/references/Celegans.fasta -simulator simlord -output DaccordCelegans -split -dazzDb reproduce_manuscript_results/correction/Daccord/CelegansDb.db -corrector daccord
cd reproduce_manuscript_results/evaluation/Daccord/
rm corrected_* reference_* uncorrected_*
