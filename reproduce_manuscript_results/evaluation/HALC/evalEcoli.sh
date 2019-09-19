#!/bin/bash

cd ../../../
python -m elector -threads 9 -uncorrected ../../longReads/simEcoli -corrected ../../correction/HALC/HALCEcoli/HALC.split.fa -reference ../../references/Ecoli.fasta -simulator simlord -output HALCEcoli -split -corrector halc
cd reproduce_manuscript_results/evaluation/HALC/
rm corrected_* reference_* uncorrected_*
