Benchmark for hybrid and self long read correction
=================================================

# Description
* This program enables the testing of different hybrid and non hybrid long read correction tools.
* It embeds a modified version of [poaV2](https://sourceforge.net/projects/poamsa/)
* Authors: Camille Marchet, Pierre Morisse, Lolita Lecompte, Antoine Limasset


# Requirements
* gcc
* Python3

# Usage

## Installation

	./install.sh

Binaries are then in ./bin

## Running the tool

	python3 benchmark.py -r reference_reads.fa -u uncorrected_reads.fa -c corrected_reads.fa

## Help

	python3 benchmark.py -h

## Toy example for testing

	python3 benchmark.py -r example/perfect_reads.fasta -c example/corrected_reads.fasta -u example/uncorrected_reads.fasta
