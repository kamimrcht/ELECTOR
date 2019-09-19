#!/bin/bash

cd ../../../
python -m elector -threads 9 -uncorrected ../../longReads/simCelegans -corrected ../../correction/HG-CoLoR/HG-CoLoRCelegans.split.fasta -reference ../../references/Celegans.fasta -simulator simlord -output HG-CoLoRCelegans -split -corrector hg-color
cd reproduce_manuscript_results/evaluation/HG-CoLoR/
rm corrected_* reference_* uncorrected_*
