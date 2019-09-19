#!/bin/bash

cd ../../../
python -m elector -threads 9 -uncorrected ../../longReads/simCelegans -corrected ../../correction/Daccord/DaccordCelegans.fasta -reference ../../references/Celegans.fasta -simulator simlord -output DaccordCelegans -split -dazzDb ../../correction/Daccord/CelegansDb.db -corrector daccord
cd reproduce_manuscript_results/evaluation/Daccord/
rm corrected_* reference_* uncorrected_*
