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
from utils import *

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
	parser.add_argument('-output', nargs='?', type=str, action="store", dest="outputDirPath", help="Name for output directory", default=None)
	parser.add_argument('-remap',  dest="remap", action='store_true', default=False, help="Perform remapping of the corrected reads to the reference")
	parser.add_argument('-assemble',  dest="assemble", action='store_true', default=False, help="Perform assembly of the corrected reads")
	parser.add_argument('-minsize', nargs='?', type=float, action="store", dest="minsize", help="Do not assess reads/fragments chose length is <= MINSIZE %% of the original read", default=10)
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
	outputDirPath = args.outputDirPath
	remap = args.remap
	assemble = args.assemble
	size_corrected_read_threshold = args.minsize / 100
	clipsNb = {}

	
	if not outputDirPath is None:
		if not os.path.exists(outputDirPath):
			os.mkdir(outputDirPath)
		else:
			printWarningMsg(outputDirPath+ " directory already exists, we will use it.")
			try:
				cmdRm = "(cd " + outputDirPath + " && rm *)"
				subprocess.check_output(['bash','-c', cmdRm])
			except subprocess.CalledProcessError:
				pass
	else:
		outputDirPath = currentDirectory
	logFile = open(outputDirPath + "/log", 'w')
	logFile.write("ELECTOR\nCommand line was:\n" + " ".join(sys.argv) + "\n")


	reportedHomopolThreshold = 5


	if perfect is not None:
		simulator = None
	if args.soft is not None:
		if args.soft == "proovread" or args.soft == "lordec" or args.soft == "nanocorr" or args.soft == "nas" or args.soft == "colormap" or args.soft == "hg-color" or args.soft == "halc" or args.soft == "pbdagcon" or args.soft == "canu" or args.soft == "lorma" or args.soft == "daccord" or args.soft == "mecat" or args.soft == "jabba" or args.soft == "fmlrc" or args.soft == "flas" or args.soft == "hercules":
			soft = args.soft
	size =  getFileReadNumber(corrected)
	if simulator is not None:
		# Récupération de la map id -> [nbLeftClips, nbRightClips] ici.
		# Si un vrai simulateur est utilisé, la map récupérée est simplement vide
		clipsNb = readAndSortFiles.processReadsForAlignment(soft, reference, uncorrected, corrected, size, split, simulator, dazzDb)
	else:
		readAndSortFiles.processReadsForAlignment(soft, perfect, uncorrected, corrected, size, split, simulator, dazzDb)

	# Check if the map is correct, ok
	# for key,val in clipsNb.items():
	# 	print (key, "=>", val[0], " ", val[1])

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
	smallReads, wronglyCorReads = alignment.getPOA(sortedCorrectedFileName, sortedRefFileName, sortedUncoFileName, args.threads, installDirectory, outputDirPath, size_corrected_read_threshold, soft)
#	alignment.getPOA(corrected, reference, uncorrected, args.threads, installDirectory, soft)
#	computeStats.outputRecallPrecision(corrected, 0, 0, soft)
	#~ print(clipsNb)
	
	nbReads, throughput, precision, recall, correctBaseRate, errorRate, smallReads, wronglyCorReads, percentGCRef, percentGCCorr, numberSplit, meanMissing, numberExtended, meanExtension, minLength, indelsubsUncorr, indelsubsCorr , truncated, ratioHomopolymer = computeStats.outputRecallPrecision(sortedCorrectedFileName, outputDirPath, logFile, smallReads, wronglyCorReads, reportedHomopolThreshold, size_corrected_read_threshold, readSizeDistribution, clipsNb, 0, 0, soft)

	#~ if simulator == "nanosim":
		#~ computeStats.outputReadSizeDistribution(uncorrected + "_reads.fasta", sortedCorrectedFileName, readSizeDistribution, outputDirPath)
	#~ elif simulator == "simlord":
		#~ computeStats.outputReadSizeDistribution(uncorrected + ".fasta", sortedCorrectedFileName, readSizeDistribution, outputDirPath)
	#~ else:
		#~ computeStats.outputReadSizeDistribution(uncorrected, sortedCorrectedFileName, readSizeDistribution, outputDirPath)
	#~ if int(numberSplit) != 0:
		#~ computeStats.outputReadSizeDistribution(sortedCorrectedFileName, readSizeDistribution, outputDirPath)

	avId=0
	cov=0
	nbContigs=0
	nbAlContig=0
	nbBreakpoints=0
	NG50=0
	NG75=0

	if remap:
		print("********** REMAPPING **********")
		logFile.write("********** REMAPPING **********\n")
		avId, cov = remappingStats.generateResults(corrected, reference, args.threads, logFile)
		print("*******************************\n")
	if assemble:
		print("********** ASSEMBLY **********")
		logFile.write("********** ASSEMBLY **********\n")
		nbContigs, nbAlContig, nbBreakpoints, NG50, NG75 = assemblyStats.generateResults(corrected, reference, args.threads, logFile)
		print("******************************")
	#TODO: inclure les nouvelles métriques dans le plot
	plotResults.generateResults(outputDirPath, installDirectory, soft, nbReads, throughput, recall, precision, correctBaseRate, errorRate, numberSplit, meanMissing, numberExtended, meanExtension, percentGCRef, percentGCCorr, smallReads, wronglyCorReads, minLength, indelsubsUncorr, indelsubsCorr, avId, cov, nbContigs, nbAlContig, nbBreakpoints, NG50, NG75 , remap, assemble, ratioHomopolymer)


if __name__ == '__main__':
	main()

