#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected reproduce_manuscript_results/longReads/simEcoli -corrected reproduce_manuscript_results/correction/MECAT/MECATEcoli.fasta -reference reproduce_manuscript_results/references/Ecoli.fasta -simulator simlord -output MECATEcoli -split -corrector mecat
cd reproduce_manuscript_results/evaluation/MECAT/
rm corrected_* reference_* uncorrected_*
