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


# sort read file by increasing order of headers
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


def lordec(threads):
	outFile = "corrected_sorted_by_lordec.fa"
	splitFile = "trimmed_split_sorted_by_lordec.fa"
	corrected_tmp = "corrected_by_lordec.fa"
	trimmedSplit_tmp = "trimmed_split_by_lordec.fa"
	logFile = open("lordec.log", 'w')
	cmdLordec = "lordec-correct -i simulatedReads.fa -2 simulatedReads_short.fa -o " + corrected_tmp + " -k 21 -s 2 -T " + str(threads)
	p = subprocessLauncher(cmdLordec, logFile, logFile)
	cmdLordecTrimSplit = "lordec-trim-split -i " + corrected_tmp + " -o " + trimmedSplit_tmp
	cmdFormatHeader = "sed -i 's/_.*//' " + trimmedSplit_tmp
	p = subprocessLauncher(cmdLordecTrimSplit, logFile, logFile)
	subprocess.check_output(['bash','-c', cmdFormatHeader])
	cmdRm = "rm *.h5"
	subprocess.check_output(['bash','-c', cmdRm])
	readAndSortFasta(corrected_tmp, outFile)
	occurenceEachRead = readAndSortFasta(trimmedSplit_tmp, splitFile)
	return (outFile, splitFile, occurenceEachRead)

#todo: formatting the headers before the sort
def colormap(threads):
	outFile = "corrected_sorted_by_colormap.fa"
	corrected_tmp = "corrected_by_colormap.fa"
	logFile = open("colormap.log", 'w')
	cmdColorMap = "runCorr.sh simulatedReads.fa simulatedReads_short.fa colorMapDir pre " + str(threads)
	p = subprocessLauncher(cmdColorMap, logFile, logFile)
	cmdColorMapOEA = "runOEA.sh colorMapDir/pre_sp.fasta simulatedReads_short.fa colorMapDir pre " + str(threads) # second pass to enhance the correction
	p = subprocessLauncher(cmdColorMapOEA, logFile, logFile)
	cmdCp = "cp colorMapDir/pre_oea.fasta " + corrected_tmp
	subprocess.check_output(['bash','-c', cmdCp])
	occurenceEachRead =  readAndSortFasta(corrected_tmp, outFile)
	return (outFile, "")

#todo: formatting the headers before the sort
def lorma(threads):
	outFile = "corrected_sorted_by_lorma.fa"
	splitFile = "trimmed_split_sorted_by_lorma.fa"
	corrected_tmp = "corrected_by_lorma.fa"
	trimmedSplit_tmp = "trimmed_split_by_lorma.fa"
	cmdcp = "cp simulatedReads.fa copy.fa"
	subprocess.check_output(['bash','-c', cmdcp])
	cmdLorma = "lorma.sh copy.fa -s -threads " + str(threads)
	p = subprocessLauncher(cmdLorma)
	occurenceEachRead = []
	try:
		cmdmv = "mv final.fasta " + corrected_tmp
		subprocess.check_output(['bash','-c', cmdmv])
		occurenceEachRead = readAndSortFasta(corrected_tmp, outFile)
	except subprocess.CalledProcessError:
		outFile = ""
	try:
		cmdmv = "mv trim.fasta " + trimmedSplit_tmp
		subprocess.check_output(['bash','-c', cmdmv])
		occurenceEachRead = readAndSortFasta(trimmedSplit_tmp, splitFile)
	except subprocess.CalledProcessError:
		splitFile = ""
	try:
		cmdmv = "mv discarded.fasta discarded_by_lorma_discarded.fa"
		subprocess.check_output(['bash','-c', cmdmv])
	except subprocess.CalledProcessError:
		pass
	
	return (outFile, splitFile, occurrenceEachRead)

#todo: provide mecat the right error model (-x)
# does not work for the moment
#~ def mecat():
	#~ outfile = "corrected_by_mecat.fa"
	#~ logFile = open("mecat.log", 'w')
	#~ cmdMecat = "mecat2pw -j 0 -d simulatedReads.fa -o candidates.txt -w . -t 4 -x 1"
	#~ p = subprocessLauncher(cmdMecat, logFile, logFile)
	#~ cmdMecat = "mecat2cns -i 0 -t 4 -x 1 candidates.txt simulatedReads.fa " + outfile
	#~ p = subprocessLauncher(cmdMecat, logFile, logFile)
	#~ return outfile


# computes msa with POA
def getPOA(corrected, reference, uncorrected, threads, soft=None):
	cmdPOA = "./bin/poa -corrected_reads_fasta " + corrected + " -reference_reads_fasta " + reference + " -uncorrected_reads_fasta " + uncorrected + " -threads " + str(threads)
	#~ cmdPOA = "./bin/poa -corrected_reads_fasta " + corrected + " -reference_reads_fasta " + reference + " -uncorrected_reads_fasta " + uncorrected
	subprocessLauncher(cmdPOA)
	if soft is not None:
		cmdMv = "mv default_output_msa.fasta msa_" + soft + ".fa"
	else:
		cmdMv = "mv default_output_msa.fasta msa.fa"
	subprocess.check_output(['bash','-c', cmdMv])



# compute recall and precision and writes output files
def outputRecallPrecision(beg=0, end=0, soft=None):
	if soft is not None:
		outProfile = open(soft + "_msa_profile.txt", 'w')
		precision, recall, missingSize = computeMetrics("msa_" + soft + ".fa", outProfile)
	else:
		outProfile = open("msa_profile.txt", 'w')
		precision, recall, missingSize = computeMetrics("msa.fa", outProfile)
	outProfile.write("\n***********SUMMARY***********\n")
	meanMissingSize = 0
	if len(missingSize) > 0:
		meanMissingSize = round(sum(missingSize)/len(missingSize),1)
	if soft is not None:
		outProfile.write(soft + ": Recall " + str(round(recall,5)) + " Precision " + str(round(precision,5)) + " Number of trimmed reads " + str(len(missingSize)) + " Mean missing size in trimmed reads " + str(meanMissingSize)+  "\n")
		print(soft + ": Recall ", round(recall,5), "Precision ", round(precision,5) , "Number of trimmed reads " , str(len(missingSize)), "Mean missing size in reads " , str(meanMissingSize))
		outProfile.write("Run in {0} seconds.".format(str(round(end-beg, 2))) + "\n") #runtime of the tool
		print(soft + " ran in {0} seconds.".format(str(round(end-beg, 2)))) #runtime of the tool
	else:
		outProfile.write("Recall " + str(round(recall,5)) + " Precision " + str(round(precision,5)) + " Number of trimmed reads " + str(len(missingSize)) + " Mean missing size in trimmed reads " + str(meanMissingSize)+  "\n")
		print("Recall ", round(recall,5), "Precision ", round(precision,5) , "Number of trimmed reads " , str(len(missingSize)), "Mean missing size in trimmed reads " , str(meanMissingSize))


def getMissingSize(reference, positionsToRemove):
	size = 0
	for position in range(positionsToRemove[0], positionsToRemove[1]+1):
		if reference[position] != ".":
			size += 1
	return size


# compute false positives, false negatives, true positives for a msa
def computeMetrics(fileName, outfile):
	msa = open(fileName, 'r')
	sumFP = []
	sumFN = []
	sumTP = []
	missingSize = []
	nbLines = 0
	lines = msa.readlines()
	readNo = 0
	while nbLines < len(lines) - 3:
		toW = ""
		if not ">" in lines[nbLines]:
			positionsToRemove = []
			reference = lines[nbLines].rstrip()
			nbLines += 2
			uncorrected =  lines[nbLines].rstrip()
			nbLines += 2
			corrected = lines[nbLines].rstrip()
			nbLines += 2
			stretches = findGapStretches(corrected)
			stretches = prolongateGapStretches(stretches)
			
			if len(stretches.keys()) > 0:
				if len(stretches.keys()) == 1:
					positionsToRemove = [next(iter(stretches.items()))[0], next(iter(stretches.items()))[1]]
					missing= getMissingSize(reference, positionsToRemove)
				else: #strange case, don't take read into account
					continue
			FN = 0 #code M
			FP = 0 #code !
			TP = 0 #code *
			position = 0
			for ntRef, ntUnco, ntResult in zip(reference, uncorrected, corrected):
				if len(positionsToRemove) == 0 or (position < positionsToRemove[0] or position > positionsToRemove[1]):
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
				position += 1
			if len(positionsToRemove) > 0:
				outfile.write(">read " + str(readNo) + "splitted_pos"+ str(positionsToRemove[0]) + ":" + str(positionsToRemove[1]) + "\n")
				outfile.write(toW + "\n")
				outfile.write("FN:" + str(FN) + " FP:" + str(FP) + " TP:" + str(TP)+"missing_size:" + str(missing) + "\n")
				missingSize.append(missing)
			else:
				outfile.write(">read " + str(readNo) + "\n")
				outfile.write(toW + "\n")
				outfile.write("FN:" + str(FN) + " FP:" + str(FP) + " TP:" + str(TP)+"\n")
				
			sumFN.append(FN)
			sumFP.append(FP)
			sumTP.append(TP)
			readNo += 1
		else:
			nbLines += 1
	
	precision = sum(sumTP)/(sum(sumTP)+sum(sumFP)) if sum(sumTP)+sum(sumFP) != 0 else 0
	recall = sum(sumTP)/(sum(sumTP)+sum(sumFN)) if  sum(sumTP)+sum(sumFN) != 0 else 0
	return precision, recall, missingSize


# get list of correctors from parameter file
def getCorrectors(parameterFileName):
	listCorrectors = []
	par = open(parameterFileName, 'r')
	lines = par.readlines()
	for line in lines:
		listCorrectors.append(line.rstrip())
	return listCorrectors



def getFileReadNumber(fileName):
	cmdGrep = """grep ">" -c """ + fileName
	val = subprocess.check_output(['bash','-c', cmdGrep])
	return int(val.decode('ascii'))


THRESH = 5
THRESH2 = 20

#find stretches of N
def findGapStretches(correctedSequence):
	prev = None
	countGap = 0
	positionsStretch = []
	for pos,ntResult in enumerate(correctedSequence):   # look for gaps in splitted/trimmed corrected read
		if prev == ".":
			if ntResult == "." and countGap > 0:  # gaps are dots in msa file
				countGap += 1
			if ntResult == "." and countGap == 0:
				countGap = 2
		if prev == None:
			if ntResult == ".":
				countGap += 1
		if ntResult != ".":
			if countGap > 0:
				positionsStretch.append([])
			countGap = 0
		if countGap >= THRESH:
			if len(positionsStretch) == 0:
				positionsStretch.append([pos-THRESH + 1, pos]) # start new stretch of gap with leftmost position
			else:
				if len(positionsStretch[-1]) == 0:
					positionsStretch[-1].extend((pos-THRESH + 1, pos))
				if len(positionsStretch[-1]) == 2:
					positionsStretch[-1][1] = pos # update position
		prev = ntResult
	return positionsStretch


#refine stretches of N:
def prolongateGapStretches(positionsStretch):
	newPositionStretch = []
	if len(positionsStretch) > 0:
		for s in positionsStretch:
			if len(s) > 0:
				if len(newPositionStretch) > 0:
					if newPositionStretch[-1][1] >= s[0] - THRESH:
						newPositionStretch[-1][1] = s[1]
				else:
					newPositionStretch.append([s[0], s[1]])
	#final stretches:
	stretch = dict()
	for s in newPositionStretch:
		if s[1] - s[0] > THRESH2:
			stretch[s[0]] = s[1]
	return stretch




def duplicateRefReads(reference, uncorrected, occurrenceEachRead, size):
	if occurrenceEachRead != [1]*size:
		refer = open(reference, 'r')
		refLines = refer.readlines()
		uncorr = open(uncorrected, 'r')
		uncoLines = uncorr.readlines()
		newRefName = "reference_duplicated.fa"
		newUncoName = "uncorrected_duplicated.fa"
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



def main():
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	# Manage command line arguments
	parser = argparse.ArgumentParser(description="Benchmark for quality assessment of long reads correctors.")
	# Define allowed options
	#~ threads = 2
	parser = argparse.ArgumentParser()
	parser.add_argument('-genome', nargs='?', type=str, action="store", dest="genomeRef", help="Reference genome file for simulation (sequence on one line)")
	parser.add_argument('-read_length', nargs='?', type=int, action="store", dest="readLen", help="Simulated read length", default=10000)
	parser.add_argument('-coverage', nargs='?', type=int, action="store", dest="coverage", help="Simulation coverage (example: 10 for 10X)", default=10)
	parser.add_argument('-error_rate', nargs='?', type=float, action="store", dest="errorRate", help="Error rate (example: 0.1 for 10 percent)", default = 0.1)
	parser.add_argument('-threads', nargs='?', type=int, action="store", dest="threads", help="Number of threads", default=2)
	parser.add_argument('-par', nargs='?', type=str, action="store", dest="parameterFile", help="parameters file: list of correctors, one per line", default=None)
	parser.add_argument('-c', nargs='?', type=str, action="store", dest="corrected", help="Fasta file with corrected reads (each read sequence on one line)")
	parser.add_argument('-u', nargs='?', type=str,  action="store", dest="uncorrected",  help="Fasta file with uncorrected reads (each read sequence on one line)")
	parser.add_argument('-r', nargs='?', type=str,  action="store", dest="reference",  help="Fasta file with reference read sequences (each read sequence on one line)")
	# get options for this run
	args = parser.parse_args()
	if (len(sys.argv) <= 1):
		parser.print_help()
		return 0
	
	corrected = ""
	uncorrected = ""
	reference = ""
	split = ""
	if args.corrected is None and args.uncorrected is None and args.reference is None:  # we simulate reads and correct them with tools incuded in the benchmark
		if args.parameterFile is None:
			print("\nParameter file is mandatory\n\n")
			parser.print_help()
			return 0
		# simulate data
		cmdSimul = "./bin/simulator " + args.genomeRef +  " " + str(args.readLen) + " " + str(args.coverage) + " " + str(args.errorRate) + " simulatedReads "
		uncorrected = "simulatedReads.fa"
		reference = "p.simulatedReads.fa"
		subprocessLauncher(cmdSimul)
		sizeUnco = getFileReadNumber(uncorrected)
		sizeRef = getFileReadNumber(reference)
		if sizeUnco == sizeRef :
			# convert fasta short reads in fastq
			fastqShortReads = open("simulatedReads_short.fq", 'w')
			cmdFq = "./bin/fa2fq simulatedReads_short.fa"
			subprocessLauncher(cmdFq, fastqShortReads)
			listCorrectors = getCorrectors(args.parameterFile)
			if len(listCorrectors) > 0:
				# launch correctors, we assume binaries are in PATH
				for soft in listCorrectors:
					print("Running", soft)
					notRun = False
					if soft == "lordec":
						beg = time.time()
						corrected, split, occurrenceEachRead = lordec(args.threads)
						reference, uncorrected = duplicateRefReads(reference, uncorrected, occurrenceEachRead, sizeRef)
						print(split)
						end = time.time()
					elif soft == "colormap":
						beg = time.time()
						corrected, split, occurrenceEachRead = colormap(args.threads)
						reference, uncorrected = duplicateRefReads(reference, uncorrected, occurrenceEachRead, sizeRef)
						end = time.time()
					elif soft == "lorma":
						beg = time.time()
						corrected, split, occurrenceEachRead = lorma(args.threads)
						end = time.time()
					else:
						notRun = True
					if not notRun:
						getPOA(split, reference, uncorrected, args.threads, soft)
						outputRecallPrecision(beg, end, soft)
		
	else: # else directly use data provided and skip simulation
		corrected = args.corrected
		uncorrected = args.uncorrected
		reference = args.reference
		getPOA(corrected, reference, uncorrected, args.threads)
		
		outputRecallPrecision()


if __name__ == '__main__':
	main()

