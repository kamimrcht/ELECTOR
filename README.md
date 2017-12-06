Benchmark for hybrid and self long read correction
=================================================

# Description
* This program enables the testing of different hybrid and non hybrid long read correction tools.
* It embeds a modified version of [poaV2](https://sourceforge.net/projects/poamsa/)
* Authors: Camille Marchet, Pierre Morisse, Lolita Lecompte, Antoine Limasset


# Requirements
* gcc and C++11
* Python3

# Usage

## Installation

	./install.sh

Binaries are then in ./bin

## Running the tool and simulate read files

	python3 benchmark.py -genome yourgenome.fa -read_length READLEN -coverage COV -error_rate RATE

## Running the tool when reference/corrected/uncorrected files are already present

	python3 benchmark.py -r reference_reads.fa -u uncorrected_reads.fa -c corrected_reads.fa

## Help

	python3 benchmark.py -h

# Toy tests

## Example 1: simulating files + bench

	python3 benchmark.py -genome example/example_reference.fasta -read_length 1000

## Example 2: directly provide read files and skip simulation

	python3 benchmark.py -r example/perfect_reads.fasta -c example/corrected_reads.fasta -u example/uncorrected_reads.fasta
