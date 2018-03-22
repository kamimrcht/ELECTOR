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
from Bio.Seq import Seq
import time
import argparse
import sys
import os
import shlex, subprocess
from subprocess import Popen, PIPE, STDOUT
import re
from os.path import basename


#Launches subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None, argstdin=None):
        args = shlex.split(cmd)
        p = subprocess.Popen(args, stdin = argstdin, stdout = argstdout, stderr = argstderr).communicate()
        return p

# format MECAT headers
def formatMecat(correctedReads, uncorrectedReads, formattedReads):
	fCor = open(correctedReads)
	fUnco = open(uncorrectedReads)
	fNewCor = open(formattedReads, 'w')

	i = 0
	hCor = fCor.readline().split("_")[0]
	hUnco = fUnco.readline()
	while hCor != "":
		sCor = fCor.readline()
		id = int(hCor.split(">")[1])
		while i != id:
			sUnco = fUnco.readline()
			hUnco = fUnco.readline()
			i = i + 1
		hCor = fCor.readline().split("_")[0]
		fNewCor.write(hUnco + sCor)

	fCor.close()
	fUnco.close()
	fNewCor.close()
	
# sort read file by increasing order of integer headers
def sortPBDCHeaders(infileName, outfileName):
	handle = open(infileName, "rU")
	l = SeqIO.parse(handle, "fasta")
	sortedList = [f for f in sorted(l, key=lambda x : int(x.id.split("/")[0]))]
	outfile = open(outfileName, 'w')
	for record in sortedList:
		outfile.write(">" + record.description+"\n")
		outfile.write(str(record.seq)+"\n")
	outfile.close()

# format daccord headers
def formatDaccord(correctedReads, uncorrectedReads, dazzDb, formattedReads):
	#dump database
	cmdDumpDb = "DBdump -rh " + dazzDb
	dumpedDb = open("daccord_dumpedDb", 'w')
	subprocessLauncher(cmdDumpDb, dumpedDb)
	fCor = open(correctedReads)
	dumpedDb.close()

	fDump = open("daccord_dumpedDb")
	lDump = fDump.readline()
	while lDump != "" and (lDump[0] == '+' or lDump[0] == '@'):
		lDump = fDump.readline()
	fUnco = open(uncorrectedReads)
	newCorrectedReads = open(formattedReads, 'w')
	curRead = 0
	hUnco = fUnco.readline()
	sUnco = fUnco.readline()

	hCor = fCor.readline()
	if hCor != '':
		hCor = hCor.split("/")[0].split(">")[1]
	while hCor != '':
		sCor = fCor.readline()
		while lDump != '' and lDump.split(" ")[1][:-1] != hCor:
			lDump = fDump.readline()
			lDump = fDump.readline()
			lDump = fDump.readline()
		lDump = fDump.readline()
		lDump = fDump.readline()
		readId = int(lDump.split(" ")[1])
		while curRead != readId:
			hUnco = fUnco.readline()
			sUnco = fUnco.readline()
			curRead = curRead + 1
		newCorrectedReads.write(hUnco + sCor)
		newHCor = fCor.readline()
		if newHCor != '':
			newHCor = newHCor.split("/")[0].split(">")[1]
		while newHCor == hCor:
			sCor = fCor.readline()
			newCorrectedReads.write(hUnco + sCor)
			newHCor = fCor.readline()
			if newHCor != '':
				newHCor = newHCor.split("/")[0].split(">")[1]
		hCor = newHCor
		lDump = fDump.readline()
	fCor.close()
	fUnco.close()
	fDump.close()
	newCorrectedReads.close()


# sort read file by increasing order of headers, return occurrence of each corrected read
def readAndSortFasta(infileName, outfileName):
	handle = open(infileName, "rU")
	l = SeqIO.parse(handle, "fasta")
	sortedList = [f for f in sorted(l, key=lambda x : (x.description))]
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
	outfile.close()
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
		for unco,ref in zip(uncoLines, refLines):
			if not ">" in ref:
				if header in occurrenceEachRead.keys():
					for times in range(occurrenceEachRead[header]):
						newRef.write(">" + header + "_" + str(times) + "\n")
						newRef.write(ref.rstrip() + "\n" )
						newUnco.write(">" + header + "_" + str(times) + "\n")
						newUnco.write(unco.rstrip() + "\n")
			else:
				header = ref.rstrip()[1:]
		return newRefName, newUncoName
	else:
		return reference, uncorrected

# format corrected reads headers
def formatHeader(corrector, correctedReads, uncorrectedReads, dazzDb, split):
	if corrector == "proovread":
		cmdFormatHeader = "sed 's/\(\.[0-9]*\)* SUBSTR.*$//g' " + correctedReads
		formattedReads = open("corrected_format_proovread.fa", 'w')
		subprocessLauncher(cmdFormatHeader, formattedReads)
		formattedReads.close()
	elif corrector == "lordec":
		if not split:
			#already formatted
			pass
		else:
			cmdFormatHeader = "sed 's/_[0-9]*$//g' " + correctedReads
			formattedReads = open("corrected_format_lordec.fa", 'w')
			subprocessLauncher(cmdFormatHeader, formattedReads)
			formattedReads.close()
	elif corrector == "nanocorr":
		cmdFormatHeader = "sed 's/_consensus$//g' " + correctedReads
		formattedReads = open("corrected_format_nanocorr.fa")
		subprocessLauncher(cmdFormatHeader, formattedReads)
		formattedReads.close()
	elif corrector == "nas":
		#already formatted
		pass
	elif corrector == "jabba":
		#TODO
		pass
	elif corrector == "colormap":
		cmdFormatHeader = "sed 's/ [0-9]* [0-9]*$//g' " + correctedReads
		formattedReads = open("corrected_format_colormap", 'w')
		subprocessLauncher(cmdFormatHeader, formattedReads)
		formattedReads.close()
	elif corrector == "hg-color":
		if not split:
			cmdFormatHeader = "sed 's/\(_-*[0-9]*\)\{4\}$//g' " + correctedReads
		else:
			cmdFormatHeader = "sed 's/\(_-*[0-9]*\)\{5\}$//g' " + correctedReads
		formattedReads = open("corrected_format_hg-color.fa", 'w')
		subprocessLauncher(cmdFormatHeader, formattedReads)
		formattedReads.close()
	elif corrector == "halc":
		if not split:
			#already formatted
			pass
		else:
			cmdFormatHeader = "sed 's/_[0-9]*$//g' " + correctedReads
			formattedReads = open("corrected_format_halc.fa", 'w')
			subprocessLauncher(cmdFormatHeader, formattedReads)
			formattedReads.close()
	elif corrector == "pbdagcon":
		sortPBDCHeaders(correctedReads, "tmp_sorted_pbdagcon.fa")
		formatDaccord("tmp_sorted_pbdagcon.fa", uncorrectedReads, dazzDb, "corrected_format_pbdagcon.fa")
	elif corrector == "canu":
		cmdFormatHeader = "sed 's/-id.*//g' " + correctedReads
		formattedReads = open("corrected_format_canu", 'w')
		subprocessLauncher(cmdFormatHeader, formattedReads)
		formattedReads.close()
	elif corrector == "lorma":
		cmdFormatHeader = "sed 's/_[0-9]*$//g' " + correctedReads
		formattedReads = open("corrected_format_lorma.fa", 'w')
		subprocessLauncher(cmdFormatHeader, formattedReads)
		formattedReads.close()
	elif corrector == "daccord":
		formatDaccord(correctedReads, uncorrectedReads, dazzDb, "corrected_format_daccord.fa")
	elif corrector == "mecat":
		formatMecat(correctedReads, uncorrectedReads, "corrected_format_mecat.fa")

def loadReference(fRef, simulator):
	f = open(fRef)
	refSeqs = {}

	if simulator == "nanosim" or simulator == "real":
		id = f.readline()[1:-1].strip().split(" ")[0].replace("_", "-")
	else:
		id = f.readline()[1:-1].strip().replace(" ", "-").replace("_", "-")
	while id != "":
		seq = f.readline()[:-1]
		refSeqs[id] = seq
		id = f.readline()[1:-1].split(" ")[0]
	f.close()
	return refSeqs

def generateRefReadsNanosim(simulatedReads, referenceGenome, referenceReads):
	fSeqs = loadReference(referenceGenome, "nanosim")
	f = open(simulatedReads)
	out = open(referenceReads, 'w')
	header = f.readline()[1:-1]
	while header != "":
		seq = f.readline()[:-1]
		id = header.split("_")
		refId = id[0]
		pos = int(id[1])
		strand = id[4]
		head = int(id[5])
		mid = int(id[6])
		tail = int(id[7])
		seq = fSeqs[refId][pos:+pos+mid]
		if strand == "R":
			seq = str(Seq(seq).reverse_complement())
		out.write(">" + header + "\n" + seq + "\n")
		header = f.readline()[1:-1]
	f.close()
	out.close()

def generateRefReadsSimLord(simulatedReads, referenceGenome, referenceReads):
	fSeqs = loadReference(referenceGenome, "simlord")
	f = open(simulatedReads)
	out = open(referenceReads, 'w')
	line = f.readline()
	while line != '' and line[0] == "@":
		line = f.readline()

	line = line.split("\t")
	while line != ['']:
		header = line[0]
		strand = int(line[1])
		refId = line[2].replace("_", "-")
		pos = int(line[3]) - 1
		cigar = line[5]
		len = int(line[8])
		nbD = sum([int(i.split("D")[0]) for i in (re.findall('\d+D', cigar))])
		nbI = sum([int(i.split("I")[0]) for i in (re.findall('\d+I', cigar))])
		len = len + nbD - nbI
		seq = fSeqs[refId][pos:pos+len+1]
		if strand == 16:
			seq = str(Seq(seq).reverse_complement())
		out.write(">" + header + "\n" + seq + "\n")
		line = f.readline().split("\t")
	f.close()
	out.close()

def generateRefReadsRealData(realReads, referenceGenome, referenceReads):
	reFile = (os.path.splitext(realReads)[0])
	cmdAl = "./minimap2/minimap2 -a " + referenceGenome + " " + realReads
	outErr = open("/dev/null", 'w')
	alFile = reFile + ".sam"
	outAl = open(alFile, 'w')
	subprocessLauncher(cmdAl, outAl, outErr)
	outAl.close()
	outErr.close()
	
	fSeqs = loadReference(referenceGenome, "real")
	f = open(alFile)
	out = open(referenceReads, 'w')
	line = f.readline()
	while line != '' and line[0] == "@":
		line = f.readline()

	line = line.split("\t")
	while line != ['']:
		if line[1] == "0" or line[1] == "16":
			header = line[0]
			strand = int(line[1])
			refId = line[2].replace("_", "-")
			pos = int(line[3]) - 1
			cigar = line[5]
			length = len(line[9])
			nbD = sum([int(i.split("D")[0]) for i in (re.findall('\d+D', cigar))])
			nbI = sum([int(i.split("I")[0]) for i in (re.findall('\d+I', cigar))])
			length = length + nbD - nbI
			seq = fSeqs[refId][pos:pos+length+1]
			if strand == 16:
				seq = str(Seq(seq).reverse_complement())
			out.write(">" + header + "\n" + seq + "\n")
		line = f.readline().split("\t")
	f.close()
	out.close()	

#Generates reference reads file (only supported for nanosim and simlord)
def convertSimulationOutputToRefFile(simulatedPrefix, referenceGenome, simulator):
	if simulator == "nanosim":
		generateRefReadsNanosim(simulatedPrefix + "_reads.fasta", referenceGenome, simulatedPrefix + "_reference.fasta")
	elif simulator == "simlord":
		cmdConv = "./bin/fq2fa " + simulatedPrefix + ".fastq"
		outFa = open(simulatedPrefix + ".fasta", 'w')
		subprocessLauncher(cmdConv, outFa)
		outFa.close()
		generateRefReadsSimLord(simulatedPrefix + ".fastq.sam", referenceGenome, simulatedPrefix + "_reference.fasta")
	else:
		generateRefReadsRealData(simulatedPrefix, referenceGenome, simulatedPrefix + "_reference.fasta")

# main function
def processReadsForAlignment(corrector, reference, uncorrected, corrected, size, split, simulator, dazzDb):
	#0- generate reference reads, if needed
	if simulator is not None:
		convertSimulationOutputToRefFile(uncorrected, reference, simulator)
	#1- correctly format the headers to be able to identify and sort the corrected reads
	if simulator == "nanosim":
		formatHeader(corrector, corrected, uncorrected + "_reads.fasta", dazzDb, split)
	elif simulator == "simlord":
		formatHeader(corrector, corrected, uncorrected + ".fasta", dazzDb, split)
	else:
		formatHeader(corrector, corrected, uncorrected, dazzDb, split)
	#2- count occurences of each corrected reads(in case of trimmed/split) and sort them
	if corrector is not None and corrector != "nas" and ((corrector != "lordec" and corrector != "halc") or split):
		newCorrectedFileName = "corrected_format_" + corrector + ".fa"
		sortedCorrectedFileName = "corrected_sorted_by_" + corrector + ".fa"
		sortedUncoFileName = "uncorrected_sorted_" + corrector + ".fa"
		newUncoFileName =  "uncorrected_sorted_duplicated_" + corrector + ".fa"
		sortedRefFileName = "reference_sorted_" + corrector + ".fa"
		newRefFileName =  "reference_sorted_duplicated_" + corrector + ".fa"
	else:
		newCorrectedFileName = corrected
		sortedCorrectedFileName = "corrected_sorted.fa"
		sortedUncoFileName = "uncorrected_sorted.fa"
		newUncoFileName =  "uncorrected_sorted_duplicated.fa"
		sortedRefFileName = "reference_sorted.fa"
		newRefFileName =  "reference_sorted_duplicated.fa"
	if simulator == "nanosim":
		readAndSortFasta(uncorrected + "_reads.fasta", sortedUncoFileName)
	elif simulator == "simlord":
		readAndSortFasta(uncorrected + ".fasta", sortedUncoFileName)
	else:
		readAndSortFasta(uncorrected, sortedUncoFileName)
	if simulator is not None:
		readAndSortFasta(uncorrected + "_reference.fasta", sortedRefFileName)
	else:
		readAndSortFasta(reference, sortedRefFileName)
	occurrenceEachRead = readAndSortFasta(newCorrectedFileName, sortedCorrectedFileName)
	#3- duplicate reference and uncorrected reads files to prepare for POA (we want as many triplets as there are corrected reads)
	duplicateRefReads(sortedRefFileName, sortedUncoFileName, occurrenceEachRead, size, newUncoFileName, newRefFileName)
