#!/bin/bash

cd ../../../
python -m elector -threads 9 -uncorrected ../../longReads/simCelegans -corrected ../../correction/MECAT/MECATCelegans.fasta -reference ../../references/Celegans.fasta -simulator simlord -output MECATCelegans -split -corrector mecat
cd reproduce_manuscript_results/evaluation/MECAT/
rm corrected_* reference_* uncorrected_*
