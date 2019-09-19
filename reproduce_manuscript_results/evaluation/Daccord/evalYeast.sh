#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected reproduce_manuscript_results/longReads/simYeast -corrected reproduce_manuscript_results/correction/Daccord/DaccordYeast.fasta -reference reproduce_manuscript_results/references/Yeast.fasta -simulator simlord -output reproduce_manuscript_results/evaluation/Daccord/DaccordYeast -split -dazzDb reproduce_manuscript_results/correction/Daccord/YeastDb.db -corrector daccord
rm corrected_* reference_* uncorrected_*
cd reproduce_manuscript_results/evaluation/Daccord/
