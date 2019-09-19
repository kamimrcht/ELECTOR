#!/bin/bash

cd ../../../
python3 -m elector -threads 16 -uncorrected reproduce_manuscript_results/longReads/simCelegans -corrected reproduce_manuscript_results/correction/Canu/CanuCelegans/CanuCelegans.correctedReads.fasta -reference reproduce_manuscript_results/references/Celegans.fasta -simulator simlord -output reproduce_manuscript_results/evaluation/Canu/CanuCelegans -corrector canu
rm corrected_* reference_* uncorrected_*
cd reproduce_manuscript_results/evaluation/Canu/
