#!/bin/bash

cd ../../../
python -m elector -threads 16 -uncorrected ../../longReads/simCelegans -corrected ../../correction/Canu/CanuCelegans/CanuCelegans.correctedReads.fasta -reference ../../references/Celegans.fasta -simulator simlord -output CanuCelegans -corrector canu
cd reproduce_manuscript_results/evaluation/Canu/
rm corrected_* reference_* uncorrected_*
