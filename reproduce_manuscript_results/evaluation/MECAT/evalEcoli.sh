#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected ../../longReads/simEcoli -corrected ../../correction/MECAT/MECATEcoli.fasta -reference ../../references/Ecoli.fasta -simulator simlord -output MECATEcoli -split -corrector mecat
cd reproduce_manuscript_results/evaluation/MECAT/
rm corrected_* reference_* uncorrected_*
