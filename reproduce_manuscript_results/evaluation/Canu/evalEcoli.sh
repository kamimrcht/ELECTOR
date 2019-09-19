#!/bin/bash

cd ../../../
python3 -m elector -threads 16 -uncorrected reproduce_manuscript_results/longReads/simEcoli -corrected reproduce_manuscript_results/correction/Canu/CanuEcoli/CanuEcoli.correctedReads.fasta -reference reproduce_manuscript_results/references/Ecoli.fasta -simulator simlord -output reproduce_manuscript_results/evaluation/Canu/CanuEcoli -corrector canu
rm corrected_* reference_* uncorrected_*
cd reproduce_manuscript_results/evaluation/Canu/
