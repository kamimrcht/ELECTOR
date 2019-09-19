#!/bin/bash

cd ../../../
python3 -m elector -threads 9 -uncorrected ../../longReads/simYeast -corrected ../../correction/MECAT/MECATYeast.fasta -reference ../../references/Yeast.fasta -simulator simlord -output MECATYeast -split -corrector mecat
cd reproduce_manuscript_results/evaluation/MECAT/
rm corrected_* reference_* uncorrected_*
