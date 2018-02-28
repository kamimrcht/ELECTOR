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
import copy

SIZE_CORRECTED_READ_THRESHOLD = 0.1

THRESH = 5
THRESH2 = 20


# store corrected reads in a list of tuples [(header, sequence)]
def getCorrectedReads(correctedReadsFileName):
	#TOVERIFY
	handle = open(correctedReadsFileName, "rU")
	l = SeqIO.parse(handle, "fasta")
	fastaTuple = {}
	for record in l:
		#fastaTuple.append((record.description, str(record.seq)))
		fastaTuple[record.description] = str(record.seq)
	return fastaTuple


# find long stretches of "." = trimmed or split reads, and return the coordinates of these regions for the corrected line of a msa
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
	stretch = dict()
	for s in positionsStretch:
		if len(s) > 0:
			if s[1] - s[0] > THRESH2:
				stretch[s[0]] = s[1]
	return stretch
	#~ return positionsStretch




# compute recall and precision and writes output files
# TODO remove runtime
def outputRecallPrecision( correctedFileName, beg=0, end=0, soft=None):
	if soft is not None:
		outProfile = open(soft + "_msa_profile.txt", 'w')
		outMetrics = open(soft + "_per_read_metrics.txt", 'w')
		precision, recall, missingSize, smallReadNumber = computeMetrics("msa_" + soft + ".fa", outProfile, outMetrics, correctedFileName)
	else:
		outProfile = open("msa_profile.txt", 'w')
		outMetrics = open("per_read_metrics.txt", 'w')
		precision, recall, missingSize, smallReadNumber = computeMetrics("msa.fa", outProfile, outMetrics, correctedFileName)
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
	print("Number of corrected reads which length is <", SIZE_CORRECTED_READ_THRESHOLD*100,"% of the original read:", smallReadNumber)
	outProfile.close()
	outMetrics.close()



def getLen(sequenceMsa):
	return len(sequenceMsa) - sequenceMsa.count('.')


# main function, compute false positives, false negatives, true positives for a msa
def computeMetrics(fileName, outMSAProfile, outPerReadMetrics, correctedFileName):
	msa = open(fileName, 'r')
	sumFP = []
	sumFN = []
	sumTP = []
	missingSize = []
	nbLines = 0
	lines = msa.readlines()
	readNo = 0
	headerNo = 0
	#TOVERIFY
#	correctedReadsList = getCorrectedReads("corrected_sorted.fa")
	correctedReadsList = getCorrectedReads(correctedFileName)
	upperCasePositions = getUpperCasePositions(correctedReadsList, lines)
	smallReadNumber = 0
	recall = 0
	precision = 0
	outPerReadMetrics.write("score metrics\n")
	while nbLines < len(lines) - 3:
		toW = ""
		if not ">" in lines[nbLines]:
			reference = lines[nbLines].rstrip()
			nbLines += 2
			uncorrected =  lines[nbLines].rstrip()
			nbLines += 2
			corrected = lines[nbLines].rstrip()
			nbLines += 2
			lenCorrected = getLen(corrected)
			lenReference = getLen(reference)
			print("************************", lenCorrected*1.0/lenReference)
			if lenCorrected*1.0/lenReference >= SIZE_CORRECTED_READ_THRESHOLD:
				stretches = findGapStretches(corrected)
				correctedPositions, missing, positionsToRemove = getCorrectedPositions(stretches, len(corrected), readNo, upperCasePositions, reference)
				FN = 0 #code M
				FP = 0 #code !
				TP = 0 #code *
				position = 0
				intervalInPositionToRemove = 0
				for ntRef, ntUnco, ntResult in zip(reference, uncorrected, corrected):
					if correctedPositions[position]:
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
				toWRead = ">read " + str(readNo)
				if len(positionsToRemove) > 0:
					for interv in positionsToRemove:
						toWRead += " splitted_pos"+ str(interv[0]) + ":" + str(interv[1]) 
					outMSAProfile.write(toWRead + "\n")
					outMSAProfile.write(toW + "\n")
					outMSAProfile.write("FN:" + str(FN) + " FP:" + str(FP) + " TP:" + str(TP)+" missing_size:" + str(missing) + "\n")
					missingSize.append(missing)
				else:
					outMSAProfile.write(toWRead + "\n")
					outMSAProfile.write(toW + "\n")
					outMSAProfile.write("FN:" + str(FN) + " FP:" + str(FP) + " TP:" + str(TP)+"\n")
					
				sumFN.append(FN)
				sumFP.append(FP)
				sumTP.append(TP)
				rec = TP / (TP + FN) if TP + FN != 0 else 0
				outPerReadMetrics.write(str(rec) + " recall\n")
				prec = TP / (TP + FP) if TP + FP != 0 else 0
				outPerReadMetrics.write(str(prec) + " precision\n")
				recall = recall + rec
				precision = precision + prec
				#outPerReadMetrics.write(str(perBaseErrorRate) + " correct_rate")
			else:
				smallReadNumber += 1
			readNo += 1
		else:
			headerNo = lines[nbLines].split(">")[1].split(" ")[0]
			nbLines += 1
	recall = recall / readNo if readNo != 0 else 0
	precision = precision / readNo if readNo != 0 else 0
	return (precision, recall, missingSize, smallReadNumber)


#TODO
# compute length missed in trimmed/split reads
def getMissingSize(reference, positionsToRemove):
	size = 0
#	print(positionsToRemove)
#	for position in range(positionsToRemove[0], positionsToRemove[1]+1):
#			if reference[position] != ".":
#				size += 1
	return size


# get the position of nt in uppercase to compute recall and precision only at these positions
def getUpperCasePositions(correctedReadsList, lines):
	upperCasePositions = [] # positions to take into account in the msa
	nbLines = 5 # starting at the first corrected sequence
	headerNo = lines[0].split(">")[1].split(" ")[0]
	while nbLines < len(lines):
		headerNo = lines[nbLines - 1].split(">")[1].split(" ")[0]
		correctedMsa = lines[nbLines].rstrip()
		#TOVERIFY
		#correctedReadSequence = correctedReadsList[int(headerNo)][1]
		correctedReadSequence = correctedReadsList[headerNo]
		posiNt = 0
		posiNtSeq = 0
		inUpper = False
		upperCasePositions.append([])
		while posiNt < len(correctedMsa):
			if posiNtSeq >= len(correctedReadSequence):
				upperCasePositions[-1].append(False)
				posiNt += 1
			else:
				nt = correctedMsa[posiNt]
				ntSeq = correctedReadSequence[posiNtSeq]
				if not ntSeq.islower():
					if nt != ".":
						inUpper = True
				else:
					inUpper = False
				
						
				if inUpper:
					upperCasePositions[-1].append(True)
				else:
					upperCasePositions[-1].append(False)
				if nt != ".":
					posiNtSeq += 1
				posiNt += 1
		nbLines += 6
		
	return upperCasePositions


# add to the uppercase positions the positions where there is no stretch of "." , i.e. all positions where recall and precision are actually computed
def getCorrectedPositions(stretches, msaLineLen, readNo, upperCasePositions, reference):
	correctedPositions = copy.copy(upperCasePositions[readNo])
	positionsToRemove = []
	missing = 0
	if len(stretches.keys()) > 0:
		if len(stretches.keys()) == 1:
			positionsToRemove = [[next(iter(stretches.items()))[0], next(iter(stretches.items()))[1]]]
			missing= getMissingSize(reference, positionsToRemove)
		else:
			for pos in stretches.keys():
				positionsToRemove.append([pos, stretches[pos]])
				missing += getMissingSize(reference, positionsToRemove[-1])
		interval = 0
		l = positionsToRemove[interval][0]
		while l  < msaLineLen:
			if l > positionsToRemove[interval][1]:
				interval += 1
				if interval < len(positionsToRemove):
					l = positionsToRemove[interval][1]
				else:
					break
			else:
				correctedPositions[l] = False
			l += 1
	return (correctedPositions, missing, positionsToRemove)
