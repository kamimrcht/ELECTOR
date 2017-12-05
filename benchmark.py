#!/usr/bin/env python3

"""*****************************************************************************
 *   Authors: Camille Marchet Antoine Limasset Pierre Morisse Lolita Lecompte
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


import time
import argparse
import sys
import os
import shlex, subprocess
from subprocess import Popen, PIPE, STDOUT
import re
import poagraph
import seqgraphalignment

# read fasta
def readfasta(infile):
    labels = []
    sequences = []

    curlabel = None
    cursequence = ""

    def updatelists():
        if len(cursequence) is not 0:
            sequences.append(cursequence)
            if curlabel is not None:
                labels.append(curlabel)
            else:
                labels.append('seq'+str(len(sequences)))

    for line in infile:
        if line[0] == ">":
            updatelists()
            cursequence = ""
            curlabel = line[1:].strip()
        else:
            cursequence += line.strip()

    updatelists()
    return list(zip(labels, sequences))


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



def computeMetrics(fastaReference, fastaUncorrected, fastaCorrected):
	sumFP = []
	sumFN = []
	sumTP = []
	for uncorrectedRead, correctedRead, referenceRead in zip(fastaUncorrected, fastaCorrected, fastaReference): #hypothesis: same order
		graph = poagraph.POAGraph(referenceRead[1], referenceRead[0]+"r")
		alignment = seqgraphalignment.SeqGraphAlignment(uncorrectedRead[1], graph, globalAlign=True)
		graph.incorporateSeqAlignment(alignment, uncorrectedRead[1], uncorrectedRead[0]+"u")
		alignment = seqgraphalignment.SeqGraphAlignment(correctedRead[1], graph, globalAlign=True)
		graph.incorporateSeqAlignment(alignment, correctedRead[1], correctedRead[0]+"c")
		alignments, consensus = graph.generateAlignmentStrings()
		toW = ""
		FN = 0 #code M
		FP = 0 #code !
		TP = 0 #code *
		for ntRef, ntUnco, ntResult in zip(alignments[0], alignments[1], alignments[2]):
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
	
	precision = sum(sumTP)/(sum(sumTP)+sum(sumFP)) if sum(sumTP)+sum(sumFP) != 0 else 0
	recall = sum(sumTP)/(sum(sumTP)+sum(sumFN)) if  sum(sumTP)+sum(sumFN) != 0 else 0
	return precision, recall


def main():
	beg = time.time()
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	# Manage command line arguments
	parser = argparse.ArgumentParser(description="Benchmark for quality assessment of long reads correctors. Usage: python3 benchmark -r perfect_reads.fa -u uncorrected_reads.fa -c corrected_reads.fa")
	# Define allowed options
	parser = argparse.ArgumentParser()
	parser.add_argument('-c', nargs='?', type=argparse.FileType('r'), action="store", dest="corrected", help="Mandatory fasta file with corrected reads (each read sequence on one line)")
	parser.add_argument('-u', nargs='?', type=argparse.FileType('r'),  action="store", dest="uncorrected",  help="Mandatory fasta file with uncorrected reads (each read sequence on one line)")
	parser.add_argument('-r', nargs='?', type=argparse.FileType('r'),  action="store", dest="reference",  help="Mandatory fasta file with reference read sequences (each read sequence on one line)")
	# get options for this run
	args = parser.parse_args()
	# format to store fasta: ('header', 'sequence')
	fastaCorrected = readfasta(args.corrected)
	fastaUncorrected = readfasta(args.uncorrected)
	fastaReference = readfasta(args.reference)
	# gets precision and recall from MSA of 3 versions of reads
	precision, recall = computeMetrics(fastaReference, fastaUncorrected, fastaCorrected)
	print("Recall:", round(recall,2), "Precision:", round(precision,2) )
	end = time.time()
	print("Run ends in {0} seconds.".format(str(round(end-beg, 2))))


		
if __name__ == '__main__':
	main()

