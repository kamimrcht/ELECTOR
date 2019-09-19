#!/bin/bash

cd ../../../
python -m elector -threads 16 -uncorrected ../../longReads/simEcoli -corrected ../../correction/Canu/CanuEcoli/CanuEcoli.correctedReads.fasta -reference ../../references/Ecoli.fasta -simulator simlord -output CanuEcoli -corrector canu
cd reproduce_manuscript_results/evaluation/Canu/
rm corrected_* reference_* uncorrected_*
