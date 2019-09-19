#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected ../../longReads/simCelegans -corrected ../../correction/LoRDEC/LordecCelegans.split.fasta -reference ../../references/Celegans.fasta -simulator simlord -output LordecCelegans -split -corrector lordec
cd reproduce_manuscript_results/evaluation/LoRDEC/
rm corrected_* reference_* uncorrected_*
