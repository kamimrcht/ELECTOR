#!/bin/bash

cd ../../../
python -m elector -threads 9 -uncorrected ../../longReads/simEcoli -corrected ../../correction/LoRDEC/LordecEcoli.split.fasta -reference ../../references/Ecoli.fasta -simulator simlord -output LordecEcoli -split -corrector lordec
cd reproduce_manuscript_results/evaluation/LoRDEC/
rm corrected_* reference_* uncorrected_*
