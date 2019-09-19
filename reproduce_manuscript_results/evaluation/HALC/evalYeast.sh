#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected reproduce_manuscript_results/longReads/simYeast -corrected reproduce_manuscript_results/correction/HALC/HALCYeast/HALC.split.fa -reference reproduce_manuscript_results/references/Yeast.fasta -simulator simlord -output HALCYeast -split -corrector halc
cd reproduce_manuscript_results/evaluation/HALC/
rm corrected_* reference_* uncorrected_*
