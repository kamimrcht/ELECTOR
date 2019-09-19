#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected reproduce_manuscript_results/longReads/simEcoli -corrected reproduce_manuscript_results/correction/HALC/HALCEcoli/HALC.split.fa -reference reproduce_manuscript_results/references/Ecoli.fasta -simulator simlord -output HALCEcoli -split -corrector halc
cd reproduce_manuscript_results/evaluation/HALC/
rm corrected_* reference_* uncorrected_*
