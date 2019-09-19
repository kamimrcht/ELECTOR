#!/bin/bash

cd ../../../
python -m elector -threads 16 -uncorrected ../../longReads/simYeast -corrected ../../correction/Canu/CanuYeast/CanuYeast.correctedReads.fasta -reference ../../references/Yeast.fasta -simulator simlord -output CanuYeast -corrector canu
cd reproduce_manuscript_results/evaluation/Canu/
rm corrected_* reference_* uncorrected_*
