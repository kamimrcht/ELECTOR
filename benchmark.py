#!/usr/bin/env python3

"""*****************************************************************************
 *   Authors: Camille Marchet  Pierre Morisse Lolita Lecompte Antoine Limasset
 *   Contact: camille.marchet@irisa.fr, IRISA/Univ Rennes/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
 *   Source: https://github.com/Malfoy/BGREAT
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



def readAndSortFasta(infileName):
	handle = open(infileName, "rU")
	l = SeqIO.parse(handle, "fasta")
	sortedList = [f for f in sorted(l, key=lambda x : x.id)]
	return sortedList

# launch subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
	args = shlex.split(cmd)
	p = subprocess.Popen(args, stdin = argstdin, stdout = argstdout, stderr = argstderr).communicate()
	return p

def checkWrittenFiles(files):
	allFilesAreOK = True
	if not os.path.isfile(files):
		print("[ERROR] There was a problem writing \"" + files + "\".")
		allFilesAreOK = False
	if not allFilesAreOK:
		dieToFatalError("One or more files could not be written.")



def computeMetrics(fileName):
	msa = open(fileName, 'r')
	sumFP = []
	sumFN = []
	sumTP = []
	nbLines = 0
	lines = msa.readlines()
	while nbLines < len(lines) - 3:
		toW = ""
		if not ">" in lines[nbLines]:
			reference = lines[nbLines].rstrip()
			nbLines += 2
			uncorrected =  lines[nbLines].rstrip()
			nbLines += 2
			corrected = lines[nbLines].rstrip()
			nbLines += 2
			FN = 0 #code M
			FP = 0 #code !
			TP = 0 #code *
			for ntRef, ntUnco, ntResult in zip(reference, uncorrected, corrected):
				if ntRef == ntUnco == ntResult:
					toW += " "
				else:
					if ntRef == ntUnco:  #no error
						if ntUnco != ntResult: #FP
							FP += 1
							toW += "!"
						# else good nt not corrected = ok
					else: #error
						if ntRef == ntResult: #error corrected
							TP += 1
							toW += "*"
						else:
							if ntUnco == ntResult: # error not corrected
								FN += 1
								toW += "M"
							else: #new error introduced by corrector
								FP += 1
								toW += "!"
			print(toW)
			print("FN:", FN, "FP:", FP, "TP:", TP)
			sumFN.append(FN)
			sumFP.append(FP)
			sumTP.append(TP)
		else:
			nbLines += 1
	
	precision = sum(sumTP)/(sum(sumTP)+sum(sumFP)) if sum(sumTP)+sum(sumFP) != 0 else 0
	recall = sum(sumTP)/(sum(sumTP)+sum(sumFN)) if  sum(sumTP)+sum(sumFN) != 0 else 0
	return precision, recall







def main():
	beg = time.time()
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	# Manage command line arguments
	parser = argparse.ArgumentParser(description="Benchmark for quality assessment of long reads correctors.")
	# Define allowed options
	parser = argparse.ArgumentParser()
	parser.add_argument('-genome', nargs='?', type=str, action="store", dest="genomeRef", help="Reference genome file for simulation (sequence on one line)")
	parser.add_argument('-read_length', nargs='?', type=int, action="store", dest="readLen", help="Simulated read length", default=10000)
	parser.add_argument('-coverage', nargs='?', type=int, action="store", dest="coverage", help="Simulation coverage (example: 10 for 10X)", default=10)
	parser.add_argument('-error_rate', nargs='?', type=float, action="store", dest="errorRate", help="Error rate (example: 0.1 for 10 percent)", default = 0.1)
	parser.add_argument('-c', nargs='?', type=str, action="store", dest="corrected", help="Fasta file with corrected reads (each read sequence on one line)")
	parser.add_argument('-u', nargs='?', type=str,  action="store", dest="uncorrected",  help="Fasta file with uncorrected reads (each read sequence on one line)")
	parser.add_argument('-r', nargs='?', type=str,  action="store", dest="reference",  help="Fasta file with reference read sequences (each read sequence on one line)")
	# get options for this run
	args = parser.parse_args()
	if (len(sys.argv) <= 1):
		parser.print_help()
		return 0
	else:
		print(args)
	corrected = ""
	uncorrected = ""
	reference = ""
	if args.corrected is None and args.uncorrected is None and args.reference is None:
		# simulate data
		cmdSimul = "./bin/simulator " + args.genomeRef +  " " + str(args.readLen) + " " + str(args.coverage) + " " + str(args.errorRate) + " simulatedReads "
		corrected = "simulatedReads.fa"
		uncorrected = "simulatedReads.fa"
		reference = "p.simulatedReads.fa"
		subprocessLauncher(cmdSimul)
	else: # else directly use data provided and skip simulation
		corrected = args.corrected
		uncorrected = args.uncorrected
		reference = args.reference
	# launch poa graph for MSA: prerequisite = all the sequences file have the same size and sequences come in the same order
	cmdPOA = "./bin/poa -corrected_reads_fasta " + corrected + " -reference_reads_fasta " + reference + " -uncorrected_reads_fasta " + uncorrected
	subprocessLauncher(cmdPOA)
	# gets precision and recall from MSA of 3 versions of reads
	precision, recall = computeMetrics("default_output_msa.fasta")
	print("Recall:", round(recall,2), "Precision:", round(precision,2))
	end = time.time()
	print("Run ends in {0} seconds.".format(str(round(end-beg, 2))))


if __name__ == '__main__':
	main()

