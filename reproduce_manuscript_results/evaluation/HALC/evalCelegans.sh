#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected reproduce_manuscript_results/longReads/simCelegans -corrected reproduce_manuscript_results/correction/HALC/HALCCelegans/HALC.split.fa -reference reproduce_manuscript_results/references/Celegans.fasta -simulator simlord -output HALCCelegans -split -corrector halc
cd reproduce_manuscript_results/evaluation/HALC/
rm corrected_* reference_* uncorrected_*
