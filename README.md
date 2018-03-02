Benchmark for hybrid and self long read correction
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

Compatible tools:

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

## Running the benchmark

	python3 benchmark.py -reference reference_reads.fa -uncorrected uncorrected_reads.fa -corrected corrected_reads.fa -threads NBTHREADS -tool TOOLNAME

* reference_reads.fa is obtained after the simulation of reads by retrieving the original sequence of each read

* uncorrected_reads.fa is the simulated reads file used to be corrected

* corrected_reads.fa is the corrected version of reads. Reads can be trimmed or not.

* the number of threads can be precised using -threads

* a tool used for correction, in lowercase and in this list (proovread, lordec, nanocorr, nas, colormap, hg-color, halc, pbdagcon, cau, lorma, daccord, mecat) can be given. In this case the pipeline will itself retrieve the correspondance between corrected, uncorrected and reference reads. If another tool is used, the user must make sure that headers in the corrected file are similar to those in the reference and uncorrected files.

## Help

	python3 benchmark.py -h


## Output

If a corrector has been given using the -tool option, the program outputs the following:

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

## Example : directly provide read files and skip simulation

	python3 benchmark.py -r example/perfect_reads.fasta -c example/corrected_reads.fasta -u example/uncorrected_reads.fasta
