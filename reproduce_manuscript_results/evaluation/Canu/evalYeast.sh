#!/bin/bash

cd ../../../
python3 -m elector -threads 16 -uncorrected reproduce_manuscript_results/longReads/simYeast -corrected reproduce_manuscript_results/correction/Canu/CanuYeast/CanuYeast.correctedReads.fasta -reference reproduce_manuscript_results/references/Yeast.fasta -simulator simlord -output CanuYeast -corrector canu
cd reproduce_manuscript_results/evaluation/Canu/
rm corrected_* reference_* uncorrected_*
