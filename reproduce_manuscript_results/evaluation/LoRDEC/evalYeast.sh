#!/bin/bash

cd ../../../
python -m elector -threads 9 -uncorrected ../../longReads/simYeast -corrected ../../correction/LoRDEC/LordecYeast.split.fasta -reference ../../references/Yeast.fasta -simulator simlord -output LordecYeast -split -corrector lordec
cd reproduce_manuscript_results/evaluation/LoRDEC/
rm corrected_* reference_* uncorrected_*
