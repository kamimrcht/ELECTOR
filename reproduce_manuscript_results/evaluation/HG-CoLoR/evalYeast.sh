#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected ../../longReads/simYeast -corrected ../../correction/HG-CoLoR/HG-CoLoRYeast.split.fasta -reference ../../references/Yeast.fasta -simulator simlord -output HG-CoLoRYeast -split -corrector hg-color
cd reproduce_manuscript_results/evaluation/HG-CoLoR/
rm corrected_* reference_* uncorrected_*
