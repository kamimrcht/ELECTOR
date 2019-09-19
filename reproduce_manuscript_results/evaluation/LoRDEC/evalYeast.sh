#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected reproduce_manuscript_results/longReads/simYeast -corrected reproduce_manuscript_results/correction/LoRDEC/LordecYeast.split.fasta -reference reproduce_manuscript_results/references/Yeast.fasta -simulator simlord -output reproduce_manuscript_results/evaluation/LoRDEC/LordecYeast -split -corrector lordec
rm corrected_* reference_* uncorrected_*
cd reproduce_manuscript_results/evaluation/LoRDEC/
