#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected reproduce_manuscript_results/longReads/simCelegans -corrected reproduce_manuscript_results/correction/MECAT/MECATCelegans.fasta -reference reproduce_manuscript_results/references/Celegans.fasta -simulator simlord -output reproduce_manuscript_results/evaluation/MECAT/MECATCelegans -split -corrector mecat
rm corrected_* reference_* uncorrected_*
cd reproduce_manuscript_results/evaluation/MECAT/
