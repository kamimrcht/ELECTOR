#!/bin/bash

cd ../../../
python -m elector -threads 9 -uncorrected ../../longReads/simYeast -corrected ../../correction/HALC/HALCYeast/HALC.split.fa -reference ../../references/Yeast.fasta -simulator simlord -output HALCYeast -split -corrector halc
cd reproduce_manuscript_results/evaluation/HALC/
rm corrected_* reference_* uncorrected_*
