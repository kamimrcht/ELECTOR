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

# sort read file by increasing order of headers, return occurrence of each corrected read
def readAndSortFasta(infileName, outfileName):
	handle = open(infileName, "rU")
	l = SeqIO.parse(handle, "fasta")
	sortedList = [f for f in sorted(l, key=lambda x : int(x.id))]
	occurrenceEachRead = dict()
	outfile = open(outfileName, 'w')
	prevHeader = ""
	for record in sortedList:
		outfile.write(">" + record.description+"\n")
		outfile.write(str(record.seq)+"\n")
		if record.description == prevHeader:
			occurrenceEachRead[record.description] += 1
		else:
			occurrenceEachRead[record.description] = 1
			prevHeader = record.description
	return occurrenceEachRead


# duplicate reference and uncorrected sequences when there are several corrected reads with the same header
def duplicateRefReads(reference, uncorrected, occurrenceEachRead, size, newUncoName, newRefName):
	if occurrenceEachRead != [1]*size:
		refer = open(reference, 'r')
		refLines = refer.readlines()
		uncorr = open(uncorrected, 'r')
		uncoLines = uncorr.readlines()
		newUnco = open(newUncoName, 'w')
		newRef = open(newRefName, 'w')
		i = 0
		for unco,ref in zip(uncoLines, refLines):
			if not ">" in ref:
				if str(i) in occurrenceEachRead.keys():
					for times in range(occurrenceEachRead[str(i)]):
						newRef.write(header + "_" + str(times) + "\n")
						newRef.write(ref.rstrip() + "\n" )
						newUnco.write(header + "_" + str(times) + "\n")
						newUnco.write(unco.rstrip() + "\n")
				i += 1
			else:
				header = ref.rstrip()
		return newRefName, newUncoName
	else:
		return reference, uncorrected

# format corrected reads headers
def formatHeader(corrector, correctedReads):
	if corrector is None or corrector == "lorma":
		cmdFormatHeader = "sed -i 's/_.*//' " + correctedReads
		subprocess.check_output(['bash','-c', cmdFormatHeader])
	
	elif corrector == "daccord":
		pass
	elif corrector == "pbdagcon":
		pass
	elif corrector == "mecat":
		pass


def convertSimulationOutputToRefFile():
	pass
	#todo : get output of the simulator and retrieve the references sequences to compare with

# main function
def processReadsForAlignment(corrector, reference, uncorrected, corrected, size, soft):
	#1- correctly format the headers to be able to identify and sort the corrected reads
	print(corrected)
	formatHeader(corrector, corrected)
	#2- count occurences of each corrected reads(in case of trimmed/split) and sort them
	if soft is not None:
		newCorrectedFileName = "corrected_sorted_by_" + soft + ".fa"
		newUncoFileName =  "uncorrected_duplicated_" + soft + ".fa"
		newRefFileName =  "reference_duplicated_" + soft + ".fa"
	else:
		newCorrectedFileName = "corrected_sorted.fa"
		newUncoFileName =  "uncorrected_duplicated.fa"
		newRefFileName =  "reference_duplicated.fa"
	occurrenceEachRead = readAndSortFasta(corrected, newCorrectedFileName)
	#3- duplicate reference and uncorrected reads files to prepare for POA (we want as many triplets as there are corrected reads)
	duplicateRefReads(reference, uncorrected, occurrenceEachRead, size, newUncoFileName, newRefFileName)

