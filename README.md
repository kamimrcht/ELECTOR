![alt text](Images/elector.png "ELECTOR.png")

ELECTOR: EvaLuation of Error Correction Tols for lOng Reads
=================================================

# Description
* This program enables the testing of different hybrid and non hybrid long read correction tools.
* It embeds a modified version of [poaV2](https://sourceforge.net/projects/poamsa/)
* Authors: Camille Marchet, Pierre Morisse, Lolita Lecompte, Antoine Limasset


# Requirements
* gcc and C++11
* Python3 and Biopython

# Usage

## Installation

	./install.sh

Binaries are then in ./bin

## Running ELECTOR

ELECTOR can be run with:

	python3 elector.py -reference referenceGenome.fa -uncorrected simulatedReadsPrefix -corrected correctedReads.fa -threads nbThreads -corrector correctorName -simulator simulatorName

where

* referenceGenome.fa is the reference genome, with one sequence per line.

* simulatedReadsPrefix is the prefix of the uncorrected reads that were simulated.

* correctedReads.fa is the corrected version of the simulated reads, with one sequence per line. The corrected reads can be trimmed and/or split, or not.

* nbThreads is the number of threads to use.

* correctorName is the corrector that was used to correct the reads. Please see the list of compatible correctors below.

* simulatorName is the simulator that was used to simulate the long reads. Please see the list of compatible simulators below.


The reference reads can also be directly provided, with:

	python3 elector.py -perfect referenceReads.fa -uncorrected uncorrectedReads.fa -corrected correctedReads.fa -threads nbThreads -corrector correctorName
	
If the corrected long reads are split, the -split option MUST be provided to ELECTOR.

## Help

	python3 elector.py -h
	
## Compatible correctors

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

If a corrector has been provided using the -tool option, the program outputs the following:

* corrector_msa_profile.txt

* msa.fa

Else, the output files must be

* msa_profile.txt

* msa.fa

The analysis of the correction per read is printed in corrector_msa_profile.txt or msa_profile.txt.
For each read, this file comes in three parts:

* A header with the identifier of the read, for instance:
	>read 45

* The alignment profile. For correct bases untouched by the corrector, we leave a blank space. For erroneous bases corrected, we print a "*". For erroneous bases left uncorrected, we print a "M". For positions were the corrector inserted a wrong base, we output a "!". For instance:

	     *                            M   M             MM                    M            M     M      M M      MM   *                         M              M  M     *         * **
	     *        **             *       *         *                      *   *   *         ***    *            *                                    *                                 
	          *       *              M        * !                          *       *                           **         *          *                       *                   **    
	               *                *  **          *** **        * ***           *     **          MM    M  *         *   *      *       *    *                  *  *           *      
	                        **      *                 * * **  *            *         *                     *     *                 ***             *  *         *              *       
	                            *                  * *                 *                           *  **      *                         *   *        *  *               M   M          
	                                 **** *        *                              M


* The false negative (FN), false positive (FP) and true positive (TP) counts for this read.

	FN = sum(M)
	
	FP = sum(!)
	
	TP = sum(*)

The last lines report of summary of the results. They provide the recall and precision of the tools and its runtime.

	recall = TP/(TP+FN)
	precision = TP/(TP+FP)

They also provide the number of trimmed reads and the mean size of missing subsequence in trimmed regions.
For trimmed reads, recall and precision are only computed on the region that is actually corrected.
If a corrected sequence is a truncated read, missing parts will also be reported in the header
For instance:

	>read 0 splitted_pos1:66 splitted_pos865:1072

indicates that a first region is missing in the beginning of the read (65 nucleotides), and a second in the end (till the 1072th last nucleotide).

The program also recalls precision, recall and runtime for each tool in the standard output.

The multiple alignment that was used to compute these results is found in msa.fa.

# Toy tests

## Example

Generate reference reads:

	python3 elector.py -reference example/example_reference.fasta -corrected example/Simlord/correctedReads.fasta -uncorrected example/Simlord/simulatedReads -simulator simlord

Directly provide reference reads:

	python3 elector.py -perfect example/Simlord/simulatedReads_reference.fasta -corrected example/Simlord/correctedReads.fasta -uncorrected example/Simlord/simulatedReads.fasta
