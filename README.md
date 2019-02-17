![alt text](Images/elector.png "ELECTOR.png")

EvaLuation of Error Correction Tools for lOng Reads
=================================================

# Description
* This program enables the testing of different hybrid and non hybrid long read correction tools.
* It embeds a modified version of [poaV2](https://sourceforge.net/projects/poamsa/)
* Authors: Camille Marchet, Pierre Morisse, Lolita Lecompte, Antoine Limasset


# Requirements
* gcc and C++11
* Python3 and Biopython
* R  

# Usage

## Installation

	git clone --recursive https://github.com/kamimrcht/ELECTOR
	
	./install.sh

Binaries are then in ./bin. ELECTOR can be run using `elector.py` in the main directory.

## Running ELECTOR

ELECTOR can be run with:

	python3 elector.py -reference referenceGenome.fa -uncorrected simulatedReadsPrefix -corrected correctedReads.fa -threads nbThreads -corrector correctorName -simulator simulatorName -output out

where

* referenceGenome.fa is the reference genome, with one sequence per line.

* simulatedReadsPrefix is the prefix of the uncorrected reads that were simulated.

* correctedReads.fa is the corrected version of the simulated reads, with one sequence per line. The corrected reads can be trimmed and/or split, or not.

* nbThreads is the number of threads to use.

* correctorName is the corrector that was used to correct the reads. Please see the list of compatible correctors below.

* simulatorName is the simulator that was used to simulate the long reads. Please see the list of compatible simulators below.

* out is a directory where to write the output


The reference reads can also be directly provided, with:

	python3 elector.py -perfect referenceReads.fa -uncorrected uncorrectedReads.fa -corrected correctedReads.fa -threads nbThreads -corrector correctorName
	
If the corrected long reads are **split**, the -split option MUST be provided to ELECTOR.

## Help

	python3 elector.py -h
	
## Current compatible correctors

* Proovread

* LoRDEC

* Nanocorr

* NaS

* CoLoRMap

* HG-CoLoR

* HALC

* PBDAGCon

* Canu

* LoRMA

* daccord

* MECAT

If one of those tools is provided with the -corrector parameter, the pipeline will itself retrieve the correspondance between corrected, uncorrected and reference reads.
If another tool is used, please do not use the -corrector parameter, and make sure that the headers of the corrected reads are similar to those of the reference and uncorrected reads.

## Compatible simulators

* SimLord

* NanoSim

If one of those tools is provided with the -simulator parameter, the pipeline will itself generate the reference reads.


## Output


## Example
Using files from example:

	cd ELECTOR
	python elector.py -uncorrected  example/uncorrected_reads_elector.fa -perfect ~/papers/ELECTOR/perfect_reads_elector.fa -corrected ~/papers/ELECTOR/corrected_reads.fa -output out -split -corrector lordec -simulator simlord
	
Output will be written in out directory.
