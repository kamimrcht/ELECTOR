#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected ../../longReads/simEcoli -corrected ../../correction/HG-CoLoR/HG-CoLoREcoli.split.fasta -reference ../../references/Ecoli.fasta -simulator simlord -output HG-CoLoREcoli -split -corrector hg-color
cd reproduce_manuscript_results/evaluation/HG-CoLoR/
rm corrected_* reference_* uncorrected_*
