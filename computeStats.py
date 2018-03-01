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
		precision, recall, missingSize, smallReadNumber, GCRateRef, GCRateCorr = computeMetrics("msa_" + soft + ".fa", outProfile, outMetrics, correctedFileName)
	else:
		outProfile = open("msa_profile.txt", 'w')
		outMetrics = open("per_read_metrics.txt", 'w')
		precision, recall, corBasesRate, missingSize, smallReadNumber, GCRateRef, GCRateCorr = computeMetrics("msa.fa", outProfile, outMetrics, correctedFileName)
	outProfile.write("\n***********SUMMARY***********\n")
	print("*********** SUMMARY ***********")
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
		print("Recall ", round(recall,5), "\nPrecision ", round(precision,5), "\nCorrect bases rate ", round(corBasesRate,5), "\nNumber of trimmed/split reads " , str(len(missingSize)), "\nMean missing size in trimmed/split reads " , str(meanMissingSize))
	print("%GC in reference reads: " + str(GCRateRef * 100) + " %GC in corrected reads: " + str(GCRateCorr * 100))
	print("Number of corrected reads which length is <", SIZE_CORRECTED_READ_THRESHOLD*100,"% of the original read:", smallReadNumber)
	outProfile.close()
	outMetrics.close()
	print("*******************************")
	return precision, recall, corBasesRate, missingSize, smallReadNumber, GCRateRef, GCRateCorr, str(len(missingSize)) , meanMissingSize



def getLen(sequenceMsa):
	return len(sequenceMsa) - sequenceMsa.count('.')

# Compute the length distribution of uncorrected and corrected reads
def outputReadSizeDistribution(uncorrectedFileName, correctedFileName, outFileName):
	unco = open(uncorrectedFileName)
	cor = open(correctedFileName)
	out = open(outFileName, 'w')
	out.write("size type\n")
	l = unco.readline()
	while l != "":
		l = unco.readline()[:-1]
		out.write(str(len(l)) + " uncorrected\n")
		l = unco.readline()
	l = cor.readline()
	while l != "":
		l = cor.readline()[:-1]
		out.write(str(len(l)) + " corrected\n")
		l = cor.readline()
	unco.close()
	cor.close()
	out.close()


def getTPFNFP(reference, uncorrected, corrected, FP, TP, FN, corBasesRates, toW, correctedPositions, positionsToRemoveBool):
	position = 0
	corBases = 0
	GCSumRef = 0
	GCSumCorr = 0
	insU = 0
	deleU = 0
	subsU = 0
	insC = 0
	deleC = 0
	subsC = 0
	for ntRef, ntUnco, ntResult in zip(reference, uncorrected, corrected):
		if ntRef.upper() == "G" or ntRef.upper() == "C":
			GCSumRef += 1
		if ntResult.upper() ==  "G" or ntResult.upper() == "C":
			GCSumCorr += 1
		#insertion deletion substitution
		if positionsToRemoveBool[position]:
			if ntRef == ".":
				if ntUnco != ".":
					insU += 1
				if ntResult != ".":
					insC += 1
			if ntRef != ".":
				if ntUnco != "." :
					if ntUnco != ntRef:
						subsU += 1
				else:
					deleU += 1
				if ntResult != "." :
					if ntResult != ntRef:
						subsC += 1
				else:
					deleC += 1
		#FP, FN, TP
		if correctedPositions[position]:
				if ntRef == ntUnco == ntResult:
					toW += " "
					corBases += 1
				else:
					if ntRef == ntUnco:  #no error
						#TODO: if inutile?
						if ntUnco != ntResult: #FP
							FP += 1
							toW += "!"
						# else good nt not corrected = ok
						else:
							corBases +=1
					else: #error
						if ntRef == ntResult: #error corrected
							TP += 1
							toW += "*"
							corBases += 1
						else:
							if ntUnco == ntResult: # error not corrected
								FN += 1
								toW += "M"
							else: #new error introduced by corrector
								FP += 1
								toW += "!"
		else:
			corBases += 1
		position += 1

	GCRateRef = round(GCSumRef * 1.0 / getLen(reference),3)
	GCRateCorr = round(GCSumCorr * 1.0 / getLen(corrected),3)
	return (FP, TP, FN, corBases, toW, GCRateRef, GCRateCorr, insU, deleU, subsU, insC, deleC, subsC)


# main function, compute false positives, false negatives, true positives for a msa
def computeMetrics(fileName, outMSAProfile, outPerReadMetrics, correctedFileName):
	msa = open(fileName, 'r')
	#TODO inutile ?
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
	corBasesRate = 0
	outPerReadMetrics.write("score metrics\n")
	prevHeader = ""
	reference = ""
	prevRef = ""
	corBasesForARead = []
	indelsubsCorr = [0,0,0] #ins del subs
	indelsubsUncorr = [0,0,0]
	FPlistForARead = []
	TPlistForARead = []
	FNlistForARead = []
	toWPrevList = []
	sameLastHeader = False
	while nbLines < len(lines) - 3:
		toW = ""
		if not ">" in lines[nbLines]:
			prevRef = reference
			reference = lines[nbLines].rstrip()
			nbLines += 2
			uncorrected =  lines[nbLines].rstrip()
			nbLines += 2
			corrected = lines[nbLines].rstrip()
			nbLines += 1
			lenCorrected = getLen(corrected)
			lenReference = getLen(reference)
			#~ if lenCorrected*1.0/lenReference >= SIZE_CORRECTED_READ_THRESHOLD:
			stretches = findGapStretches(corrected)
			correctedPositions, positionsToRemove, positionsToRemoveBool = getCorrectedPositions(stretches, len(corrected), readNo, upperCasePositions, reference)
			#~ print("corrected positions", correctedPositions)
			corBases = 0
			FN = 0 #code M
			FP = 0 #code !
			TP = 0 #code *
			corBasesRateForARead = 0
			position = 0
			intervalInPositionToRemove = 0
			FP, TP, FN, corBases, toW, GCRateRef, GCRateCorr, insU, deleU, subsU, insC, deleC, subsC = getTPFNFP(reference, uncorrected, corrected, FP, TP, FN, corBases, toW, correctedPositions, positionsToRemoveBool)
			indelsubsCorr[0] += insC
			indelsubsCorr[1] += deleC
			indelsubsCorr[2] += subsC
			indelsubsUncorr[0] += insU
			indelsubsUncorr[1] += deleU
			indelsubsUncorr[2] += subsU
			if headerNo == prevHeader:
				sameLastHeader = True
				corBasesForARead.append(corBases)
				FPlistForARead.append(FP)
				TPlistForARead.append(TP)
				FNlistForARead.append(FN)
				toWPrevList.append(toW)
				#~ indelsubsCorr[0] += insC
				#~ indelsubsCorr[1] += deleC
				#~ indelsubsCorr[2] += subsC
				#~ indelsubsUncorr[0] += insU
				#~ indelsubsUncorr[1] += deleU
				#~ indelsubsUncorr[2] += subsU
				#todo change
				toWPrev = toWPrevList[-1]
				positionsToRemovePrev = positionsToRemove
				positionsToRemoveBoolPrev = [any(tup) for tup in zip(positionsToRemoveBoolPrev, positionsToRemoveBool)] #logical OR
				prevCorrectedPositions = [any(tup) for tup in zip(prevCorrectedPositions, correctedPositions)] #logical OR
				#~ missingPrev = positionsToRemoveBoolPrev.count(False)
				missingPrevPositions = [i for i,x in enumerate(positionsToRemoveBoolPrev) if x == False]
				missingPrev = sum([1 for x in missingPrevPositions if reference[x] != "."])
			else:
				prevCorrectedPositions = correctedPositions
				if prevCorrectedPositions.count(True)*1.0/lenReference >= SIZE_CORRECTED_READ_THRESHOLD:
					if prevHeader != "":
						corBasesPrev = sum(corBasesForARead)*1.0/len(corBasesForARead) if len(corBasesForARead) > 0 else 0
						FNPrev = sum(FNlistForARead)*1.0/len(FNlistForARead) if len(FNlistForARead) > 0 else 0
						FPPrev = sum(FPlistForARead)*1.0/len(FPlistForARead) if len(FPlistForARead) > 0 else 0
						TPPrev = sum(TPlistForARead)*1.0/len(TPlistForARead) if len(TPlistForARead) > 0 else 0
						toWReadPrev = ">read " + str(readNoPrev)
						#todo change
						#~ if len(positionsToRemovePrev) > 0:
							#~ for interv in positionsToRemovePrev:
								#~ toWReadPrev += " splitted_pos"+ str(interv[0]) + ":" + str(interv[1]) 
						if missingPrev != 0:
							missingSize.append(missingPrev)
						# append ?
						sumFN.append(FNPrev)
						sumFP.append(FPPrev)
						sumTP.append(TPPrev)
						#~ indelsubsCorr[0] = insC
						#~ indelsubsCorr[1] = deleC
						#~ indelsubsCorr[2] = subsC
						#~ indelsubsUncorr[0] = insU
						#~ indelsubsUncorr[1] = deleU
						#~ indelsubsUncorr[2] = subsU
					else: #first triplet
						#~ indelsubsCorr[0] += insC
						#~ indelsubsCorr[1] += deleC
						#~ indelsubsCorr[2] += subsC
						#~ indelsubsUncorr[0] += insU
						#~ indelsubsUncorr[1] += deleU
						#~ indelsubsUncorr[2] += subsU
						positionsToRemovePrev = positionsToRemove
						positionsToRemoveBoolPrev = positionsToRemoveBool
						prevCorrectedPositions = correctedPositions
						#~ missingPrev = positionsToRemoveBoolPrev.count(False)
						missingPrevPositions = [i for i,x in enumerate(positionsToRemoveBoolPrev) if x == False]
						missingPrev = sum([1 for x in missingPrevPositions if reference[x] != "."])
				else:
					smallReadNumber += 1
				if FPlistForARead != [] or TPlistForARead != [] or FNlistForARead != [] or corBasesForARead != []:
					corBasesSum = sum(corBasesForARead)
					FNsum = sum(FNlistForARead)
					TPsum = sum(TPlistForARead)
					FPsum = sum (FPlistForARead)
					rec = TPsum / (TPsum + FNsum) if (TPsum + FNsum) != 0 else 0
					prec = TPsum / (TPsum + FPsum) if (TPsum + FPsum) != 0 else 0
					corBRate = corBasesSum / len(prevRef)
					outPerReadMetrics.write(str(rec) + " recall\n")
					outPerReadMetrics.write(str(prec) + " precision\n")
					outPerReadMetrics.write(str(corBRate) + " correct_rate\n")
					recall += rec
					precision += prec
					corBasesRate += corBRate
				corBasesForARead = [corBases]
				FPlistForARead = [FP]
				TPlistForARead = [TP]
				FNlistForARead = [FN]
				toWPrevList = [toW]
				positionsToRemovePrev = positionsToRemove
				positionsToRemoveBoolPrev = positionsToRemoveBool
				prevCorrectedPositions = correctedPositions
			sumFN.append(FN)
			sumFP.append(FP)
			sumTP.append(TP)
#			rec = TP / (TP + FN) if TP + FN != 0 else 0
#			outPerReadMetrics.write(str(rec) + " recall\n")
#			prec = TP / (TP + FP) if TP + FP != 0 else 0
#			outPerReadMetrics.write(str(prec) + " precision\n")
#			outPerReadMetrics.write(str(corBasesRateForARead) + " correct_rate\n")
#			recall = recall + rec
#			precision = precision + prec
#			corBasesRate = corBasesRate + corBasesRateForARead
			prevHeader = headerNo
			readNoPrev = readNo
			readNo += 1
		else:
			headerNo = lines[nbLines].split(">")[1].split(" ")[0]
			nbLines += 1
	if sameLastHeader:
		if prevCorrectedPositions.count(True)*1.0/lenReference >= SIZE_CORRECTED_READ_THRESHOLD:
			FNPrev = round(sum(FNlistForARead)*1.0/len(FNlistForARead),1) if len(FNlistForARead) > 0 else 0
			FPPrev = round(sum(FPlistForARead)*1.0/len(FPlistForARead),1) if len(FPlistForARead) > 0 else 0
			TPPrev = round(sum(TPlistForARead)*1.0/len(TPlistForARead),1) if len(TPlistForARead) > 0 else 0
			toWReadPrev = ">read " + str(readNoPrev)
			#todo change
			#~ if len(positionsToRemovePrev) > 0:
				#~ for interv in positionsToRemovePrev:
					#~ toWReadPrev += " splitted_pos"+ str(interv[0]) + ":" + str(interv[1]) 
			if missingPrev != 0:
				missingSize.append(missingPrev)
			sumFN.append(FNPrev)
			sumFP.append(FPPrev)
			sumTP.append(TPPrev)
		else:
			smallReadNumber += 1

	if FPlistForARead != [] or TPlistForARead != [] or FNlistForARead != [] or corBasesForARead != []:
		corBasesSum = sum(corBasesForARead)
		FNsum = sum(FNlistForARead)
		TPsum = sum(TPlistForARead)
		FPsum = sum (FPlistForARead)
		rec = TPsum / (TPsum + FNsum) if (TPsum + FNsum) != 0 else 0
		prec = TPsum / (TPsum + FPsum) if (TPsum + FPsum) != 0 else 0
		corBRate = corBasesSum / len(reference)
		outPerReadMetrics.write(str(rec) + " recall\n")
		outPerReadMetrics.write(str(prec) + " precision\n")
		outPerReadMetrics.write(str(corBRate) + " correct_rate\n")
		recall += rec
		precision += prec
		corBasesRate += corBRate

	recall = recall / readNo if readNo != 0 else 0
	precision = precision / readNo if readNo != 0 else 0
	corBasesRate = corBasesRate / readNo if readNo != 0 else 0
	#~ print(indelsubsCorr[0],indelsubsCorr[1],indelsubsCorr[2],indelsubsUncorr[0],indelsubsUncorr[1],indelsubsUncorr[2])
	return (precision, recall, corBasesRate, missingSize, smallReadNumber, GCRateRef, GCRateCorr,indelsubsUncorr, indelsubsCorr)




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
	positionsToRemoveBool = [True] * len(correctedPositions)
	positionsToRemove = []
	#~ missing = 0
	if len(stretches.keys()) > 0:
		if len(stretches.keys()) == 1:
			positionsToRemove = [[next(iter(stretches.items()))[0], next(iter(stretches.items()))[1]]]
			#~ missing= getMissingSize(reference, positionsToRemove)
		else:
			for pos in stretches.keys():
				positionsToRemove.append([pos, stretches[pos]])
				#~ missing += getMissingSize(reference, positionsToRemove[-1])
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
	for interv in positionsToRemove:
		for i in range(interv[0], interv[1] + 1):
			positionsToRemoveBool[i] = False
	return (correctedPositions, positionsToRemove, positionsToRemoveBool)
