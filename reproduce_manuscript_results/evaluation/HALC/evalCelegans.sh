#!/bin/bash

cd ../../../
python -m elector -threads 9 -uncorrected ../../longReads/simCelegans -corrected ../../correction/HALC/HALCCelegans/HALC.split.fa -reference ../../references/Celegans.fasta -simulator simlord -output HALCCelegans -split -corrector halc
cd reproduce_manuscript_results/evaluation/HALC/
rm corrected_* reference_* uncorrected_*
