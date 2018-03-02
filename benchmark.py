#!/usr/bin/env python3

"""*****************************************************************************
 *   Authors: Camille Marchet  Pierre Morisse Antoine Limasset
 *   Contact: camille.marchet@irisa.fr, IRISA/Univ Rennes/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
 *   Source: https://github.com/kamimrcht/benchmark-long-read-correction
 *
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************"""

from Bio import SeqIO
import time
import argparse
import sys
import os
import shlex, subprocess
from subprocess import Popen, PIPE, STDOUT
import re 
#~ import readAndSortFiles
#~ import computeStats
import alignment
import computeStats
import readAndSortFiles
import plotResults
import remappingStats
import assemblyStats


# count the number of reads in a file
def getFileReadNumber(fileName):
	cmdGrep = """grep ">" -c """ + fileName
	val = subprocess.check_output(['bash','-c', cmdGrep])
	return int(val.decode('ascii'))



def main():
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	installDirectory = os.path.dirname(os.path.realpath(__file__))
	# Manage command line arguments
	parser = argparse.ArgumentParser(description="Benchmark for quality assessment of long reads correctors.")
	# Define allowed options
	corrected = ""
	uncorrected = ""
	perfect = ""
	reference = ""
	simulator = ""
	parser = argparse.ArgumentParser()
	parser.add_argument('-threads', nargs='?', type=int, action="store", dest="threads", help="Number of threads", default=2)
	parser.add_argument('-corrected', nargs='?', type=str, action="store", dest="corrected", help="Fasta file with corrected reads (each read sequence on one line)")
	parser.add_argument('-split',  dest="split", action='store_true', default=False, help="Corrected reads are split")
	parser.add_argument('-uncorrected', nargs='?', type=str,  action="store", dest="uncorrected",  help="Prefix of the reads simulation files")
	parser.add_argument('-perfect', nargs='?', type=str, action="store", dest="perfect", help="Fasta file with reference read sequences (each read sequence on one line)")
	parser.add_argument('-reference', nargs='?', type=str,  action="store", dest="reference",  help="Fasta file with reference genome sequences (each sequence on one line)")
	parser.add_argument('-simulator', nargs='?', type=str, action="store", dest="simulator", help="Tool used for the simulation of the long reads (either nanosim or simlord)")
	parser.add_argument('-corrector', nargs='?', type=str,  action="store", dest="soft",  help="Corrector used (lowercase, in this list: lorma, mecat, pbdagcon, daccord). If no corrector name is provided, make sure the read's headers are correctly formatted (i.e. they correspond to those of uncorrected and reference files)")
	parser.add_argument('-dazzDb', nargs='?', type=str, action="store", dest="dazzDb", help="Reads database used for the correction, if the reads were corrected with Daccord or PBDagCon")
	# get options for this run
	args = parser.parse_args()
	if (len(sys.argv) <= 1):
		parser.print_help()
		return 0

	corrected = args.corrected
	uncorrected = args.uncorrected
	perfect = args.perfect
	reference = args.reference
	split = args.split
	soft = None
	dazzDb = args.dazzDb
	simulator = args.simulator
	if perfect is not None:
		simulator = None
	if args.soft is not None:
		if args.soft == "proovread" or args.soft == "lordec" or args.soft == "nanocorr" or args.soft == "nas" or args.soft == "colormap" or args.soft == "hg-color" or args.soft == "halc" or args.soft == "pbdagcon" or args.soft == "canu" or args.soft == "lorma" or args.soft == "daccord" or args.soft == "mecat":
			soft = args.soft
	size =  getFileReadNumber(corrected)
	if simulator is not None:
		readAndSortFiles.processReadsForAlignment(soft, reference, uncorrected, corrected, size, split, simulator, dazzDb)
	else:
		readAndSortFiles.processReadsForAlignment(soft, perfect, uncorrected, corrected, size, split, simulator, dazzDb)
	#TOVERIFY
	if soft is not None:
		sortedCorrectedFileName = "corrected_sorted_by_" + soft + ".fa"
		sortedUncoFileName =  "uncorrected_sorted_duplicated_" + soft + ".fa"
		sortedRefFileName =  "reference_sorted_duplicated_" + soft + ".fa"
		readSizeDistribution = soft + "_read_size_distribution.txt"
	else:
		sortedCorrectedFileName = "corrected_sorted.fa"
		sortedUncoFileName =  "uncorrected_sorted_duplicated.fa"
		sortedRefFileName =  "reference_sorted_duplicated.fa"
		readSizeDistribution = "read_size_distribution.txt"
	alignment.getPOA(sortedCorrectedFileName, sortedRefFileName, sortedUncoFileName, args.threads, installDirectory, soft)
#	alignment.getPOA(corrected, reference, uncorrected, args.threads, installDirectory, soft)
#	computeStats.outputRecallPrecision(corrected, 0, 0, soft)
	precision, recall, correctBaseRate, meanMissing, smallReads, percentGCRef, percentGCCorr, numberSplit, meanMissing, indelsubsUncorr, indelsubsCorr  = computeStats.outputRecallPrecision(sortedCorrectedFileName, 0, 0, soft)

	if simulator == "nanosim":
		computeStats.outputReadSizeDistribution(uncorrected + "_reads.fasta", sortedCorrectedFileName, readSizeDistribution)
	elif simulator == "simlord":
		computeStats.outputReadSizeDistribution(uncorrected + ".fasta", sortedCorrectedFileName, readSizeDistribution)
	else:
		computeStats.outputReadSizeDistribution(uncorrected, sortedCorrectedFileName, readSizeDistribution)

	if reference is not None:
		print("********** REMAPPING **********")
		avId, cov = remappingStats.generateResults(corrected, reference, args.threads)
		print("*******************************")
		print("********** ASSEMBLY **********")
		nbContigs, nbAlContig, nbBreakpoints, NG50, NG75 = assemblyStats.generateResults(corrected, reference, args.threads)
		print("******************************")
	plotResults.generateResults(currentDirectory, installDirectory, soft, recall, precision, correctBaseRate, numberSplit, meanMissing, percentGCRef, percentGCCorr, smallReads, indelsubsUncorr, indelsubsCorr, avId, cov, nbContigs, nbAlContig, nbBreakpoints, NG50, NG75 )


if __name__ == '__main__':
	main()

