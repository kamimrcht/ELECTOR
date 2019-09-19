#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected reproduce_manuscript_results/longReads/simCelegans -corrected reproduce_manuscript_results/correction/HG-CoLoR/HG-CoLoRCelegans.split.fasta -reference reproduce_manuscript_results/references/Celegans.fasta -simulator simlord -output HG-CoLoRCelegans -split -corrector hg-color
cd reproduce_manuscript_results/evaluation/HG-CoLoR/
rm corrected_* reference_* uncorrected_*
