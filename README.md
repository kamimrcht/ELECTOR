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

Any tool you wish to test must be installed and in your PATH.

To add a tool in your PATH:

	PATH=/path/to/tool/binary:$PATH

Compatible tools:

* LoRDEC

* ColorMap

## Running the tool and simulate read files

	python3 benchmark.py -genome yourgenome.fa -read_length READLEN -coverage COV -error_rate RATE -par PARAMETERS_FILE -threads NBTHREADS

* the mandatory genome file must be in fasta format, with one line per sequences

* optional read_length, coverage and error rate modify the length, coverage and error rate of long reads (for instance -read_length 10000 for size 10000, -coverage 10 for 10X, -error_rate 0.1 for 10 percent)

* the mandatory parameter file is a text file with name of the correctors to be run, one by line, in lower case

* optionnally you can ask for more threads (default is 2), this will be used both for accelerating multiple alignment and correctors runs

The simulation produces 3 files in the directory:

* simulatedReads.fa is the long erroneous read file

* p.simulatedReads.fa is the long genomic read file (same sequences, in same order than in simulatedReads.fa, without errors)

* simulatedReads_short.fa is the short reads file

## Running the tool when reference/corrected/uncorrected files are already present

	python3 benchmark.py -r reference_reads.fa -u uncorrected_reads.fa -c corrected_reads.fa -threads NBTHREADS

## Help

	python3 benchmark.py -h


## Output

For each corrector, the program outputs two files

* corrector.log

* corrector_msa_profile.txt

The log of the corrector's run is printed in corrector.log file.
The analysis of the correction per read is printed in corrector_msa_profile.txt.
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

The program also recalls precision, recall and runtime for each tool in the standard output.

# Toy tests

## Example 1: simulating files + bench

	python3 benchmark.py -genome example/example_reference.fasta -read_length 1000 -par example/correctors.par -threads 4

## Example 2: directly provide read files and skip simulation

	python3 benchmark.py -r example/perfect_reads.fasta -c example/corrected_reads.fasta -u example/uncorrected_reads.fasta
