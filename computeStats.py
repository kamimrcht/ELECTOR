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
	fastaTuple = dict()
	#~ read = ""
	with open(correctedReadsFileName) as fileIn:
		for line in fileIn:
			if ">" in line:
				header = line.rstrip().split(' ')[0][1:]
			else:
				read = line.rstrip()
				if read != "":
					if header not in fastaTuple.keys():
						fastaTuple[header] = [read]
					else:
						fastaTuple[header].append(read)
	return fastaTuple


# find long stretches of "." = trimmed or split reads, and return the coordinates of these regions for the corrected line of a msa
# if reference is not a "." too, else it means it is a gap opened by the uncorrected part
def findGapStretches(correctedSequence, referenceSequence):
	prev = None
	countGap = 0
	countGapRef = 0
	positionsStretch = []
	pos = 0
	#~ for pos,ntResult in enumerate(correctedSequence):   # look for gaps in splitted/trimmed corrected read
	for ntRef,ntResult in zip(referenceSequence, correctedSequence):   # look for gaps in splitted/trimmed corrected read
		if prev == ".":
			if ntResult == "." and countGap > 0:  # gaps are dots in msa file
				countGap += 1
			if ntResult == "." and countGap == 0:
				countGap = 2
			if ntRef == "." and countGapRef > 0:  # gaps are dots in msa file
				countGapRef += 1
			if ntRef == "." and countGapRef == 0:
				countGapRef = 2
		if prev == None:
			if ntResult == ".":
				countGap += 1
			if ntRef == ".":
				countGapRef += 1
		if ntResult != ".":
			if countGap > 0:
				positionsStretch.append([])
			countGap = 0
		if ntRef != ".":
			countGapRef = 0
		if countGap >= THRESH:
			if countGapRef < THRESH: #if countGapRef>=THRESH it means that a gap is openened both in ref and corrected because of the uncorrected seq in the msa, so this is not a trimmed zone
				if len(positionsStretch) == 0:
					positionsStretch.append([pos-THRESH + 1, pos]) # start new stretch of gap with leftmost position
				else:
					if len(positionsStretch[-1]) == 0:
						positionsStretch[-1].extend((pos-THRESH + 1, pos))
					if len(positionsStretch[-1]) == 2:
						positionsStretch[-1][1] = pos # update position
		prev = ntResult
		pos += 1
	stretch = dict()
	for s in positionsStretch:
		if len(s) > 0:
			if s[1] - s[0] > THRESH2:
				stretch[s[0]] = s[1]
	return stretch
	#~ return positionsStretch




# compute recall and precision and writes output files
# TODO remove runtime
def outputRecallPrecision( correctedFileName, outDir, logFile, smallReadNumber, beg=0, end=0,  soft=None):
	if soft is not None:
		outProfile = open(outDir + "/" + soft + "_msa_profile.txt", 'w')
		outMetrics = open(outDir + "/" + soft + "_per_read_metrics.txt", 'w')
		precision, recall, corBasesRate, missingSize,  GCRateRef, GCRateCorr,  indelsubsUncorr, indelsubsCorr = computeMetrics(outDir + "/msa_" + soft + ".fa", outProfile, outMetrics, correctedFileName)
	else:
		outProfile = open(outDir + "/msa_profile.txt", 'w')
		outMetrics = open(outDir + "/per_read_metrics.txt", 'w')
		precision, recall, corBasesRate, missingSize,  GCRateRef, GCRateCorr, indelsubsUncorr, indelsubsCorr = computeMetrics(outDir + "/msa.fa", outProfile, outMetrics, correctedFileName)
	outProfile.write("\n***********SUMMARY***********\n")
	print("*********** SUMMARY ***********")
	meanMissingSize = 0
	if len(missingSize) > 0:
		meanMissingSize = round(sum(missingSize)/len(missingSize),1)
	if soft is not None:
		outProfile.write(soft + "Recall " + str(round(recall,5)) + " Precision " + str(round(precision,5)) + " Number of trimmed reads " + str(len(missingSize)) + " Mean missing size in trimmed reads " + str(meanMissingSize)+  "\n")
		print(soft + "\nRecall :", round(recall,5), "\nPrecision :", round(precision,5), "\nCorrect bases rate :", round(corBasesRate,5), "\nNumber of trimmed/split reads :" , str(len(missingSize)), "\nMean missing size in trimmed/split reads :" , str(meanMissingSize))
		outProfile.write("Run in {0} seconds.".format(str(round(end-beg, 2))) + "\n") #runtime of the tool
		print(soft + " ran in {0} seconds.".format(str(round(end-beg, 2)))) #runtime of the tool
	else:
		outProfile.write("Recall " + str(round(recall,5)) + " Precision " + str(round(precision,5)) + " Number of trimmed reads " + str(len(missingSize)) + " Mean missing size in trimmed reads " + str(meanMissingSize)+  "\n")
		print("Recall :", round(recall,5), "\nPrecision :", round(precision,5), "\nCorrect bases rate :", round(corBasesRate,5), "\nNumber of trimmed/split reads :" , str(len(missingSize)), "\nMean missing size in trimmed/split reads :" , str(meanMissingSize))
	print("%GC in reference reads : " + str(GCRateRef * 100) + "\n%GC in corrected reads : " + str(GCRateCorr * 100))
	print("Number of corrected reads which length is <", SIZE_CORRECTED_READ_THRESHOLD*100,"% of the original read :", smallReadNumber)
	outProfile.close()
	outMetrics.close()
	print("Number of insertions in uncorrected : " + str(indelsubsUncorr[0]) +"\nNumber of insertions in corrected : " + str(indelsubsCorr[0]) )
	print("Number of deletions in uncorrected : " + str(indelsubsUncorr[1]) +"\nNumber of deletions in corrected : " + str(indelsubsCorr[1]) )
	print("Number of substitutions in uncorrected : " + str(indelsubsUncorr[2]) +"\nNumber of substitutions in corrected : " + str(indelsubsCorr[2]) )
	
	logFile.write("*********** SUMMARY ***********\nRecall :" + str(round(recall,5)) + "\nPrecision :" + str(round(precision,5)) + "\nCorrect bases rate :" + str(round(corBasesRate,5)) + "\nNumber of trimmed/split reads :" + str(len(missingSize)) + "\nMean missing size in trimmed/split reads :" + str(meanMissingSize) + "\n%GC in reference reads : " + str(GCRateRef * 100) + "\n%GC in corrected reads : " + str(GCRateCorr * 100) + "\nNumber of corrected reads which length is <" + str(SIZE_CORRECTED_READ_THRESHOLD*100) + "% of the original read :" + str(smallReadNumber) + "\nNumber of insertions in uncorrected : " + str(indelsubsUncorr[0]) +"\nNumber of insertions in corrected : " + str(indelsubsCorr[0]) + "\nNumber of deletions in uncorrected : " + str(indelsubsUncorr[1]) +"\nNumber of deletions in corrected : " + str(indelsubsCorr[1]) + "\nNumber of substitutions in uncorrected : " + str(indelsubsUncorr[2]) +"\nNumber of substitutions in corrected : " + str(indelsubsCorr[2]) +"\n")
	return precision, recall, corBasesRate, missingSize, smallReadNumber, GCRateRef, GCRateCorr, str(len(missingSize)) , meanMissingSize,  indelsubsUncorr, indelsubsCorr



def getLen(sequenceMsa):
	return len(sequenceMsa) - sequenceMsa.count('.')

# Compute the length distribution of uncorrected and corrected reads
def outputReadSizeDistribution(uncorrectedFileName, correctedFileName, outFileName, outDir):
	unco = open(uncorrectedFileName)
	cor = open(correctedFileName)
	out = open(outDir + "/" + outFileName, 'w')
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


# compute ins, del, subs
def indels(ntRef, ntUnco, ntResult,  existingCorrectedPositions, position, insU, deleU, subsU, insC, deleC, subsC):
	#compute indels in uncorrected reads
	if existingCorrectedPositions[position]:
		if ntUnco != ntRef:
			if ntRef == ".":
				insU += 1
			else:
				if ntUnco != "." :
					subsU += 1
				else:
					deleU += 1
	#compute only indels in parts of the MSA that actually correspond to a portion that exist in the corrected read
	if existingCorrectedPositions[position]:
		if ntResult != ntRef:
			if ntRef == ".":
				insC += 1
			else:
				if ntResult != "." :
					subsC += 1
				else:
					deleC += 1
	return insU, deleU, subsU, insC, deleC, subsC

# compute fp, tp, fn
def getCorrectionAtEachPosition(ntRef, ntUnco, ntResult, correctedPositions, existingCorrectedPositions, position, corBases, uncorBase, FP, FN, TP):
	#compute FN,FP,TP only on corrected parts
	if correctedPositions[position]:
		if ntRef == ntUnco == ntResult:
			corBases += 1
		else:
			if ntRef == ntUnco:  #no error
				if ntUnco != ntResult: #FP
					FP += 1
					uncorBase += 1
				# else good nt not corrected = ok
				else:
					corBases +=1
			else: #error
				if ntRef == ntResult: #error corrected
					TP += 1
					corBases += 1
				else:
					if ntUnco == ntResult: # error not corrected
						FN += 1
					else: #new error introduced by corrector
						FP += 1
					uncorBase += 1
	# correct base rate is computed everywhere
	else:
		if existingCorrectedPositions[position] :
			if ntRef == ntResult:
				corBases += 1
			else:
				uncorBase += 1
	return corBases, uncorBase, FP, FN, TP


# get insertion deletion substitution FP, FN, TP and GC rates for a triplet
def getTPFNFP(reference, uncorrected, corrected,  correctedPositions, existingCorrectedPositions):
	position = 0
	corBases = 0
	uncorBases = 0
	FN = 0
	FP = 0
	TP = 0
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
		insU, deleU, subsU, insC, deleC, subsC = indels(ntRef, ntUnco, ntResult, existingCorrectedPositions, position, insU, deleU, subsU, insC, deleC, subsC)
		#FP, FN, TP
		corBases, uncorBases, FP, FN, TP = getCorrectionAtEachPosition(ntRef, ntUnco, ntResult, correctedPositions,  existingCorrectedPositions, position,  corBases, uncorBases, FP, FN, TP)
		position += 1
	GCRateRef = round(GCSumRef * 1.0 / getLen(reference),3)
	GCRateCorr = round(GCSumCorr * 1.0 / getLen(corrected),3)
	return FP, TP, FN, corBases, uncorBases, GCRateRef, GCRateCorr, insU, deleU, subsU, insC, deleC, subsC


def outputMetrics(recall, precision, corBasesRate, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, lenReference, existingCorrectedPositionsInThisRead, FPlistForARead, TPlistForARead, FNlistForARead, corBasesForARead, uncorBasesForARead, GCRateRefRead, GCRateCorrRead, nbReadsToDivide):
	missingInRead = existingCorrectedPositionsInThisRead.count(False) # final sum 
	if FPlistForARead != [] or TPlistForARead != [] or FNlistForARead != [] or corBasesForARead != []:
		FNsum = sum(FNlistForARead)
		TPsum = sum(TPlistForARead)
		FPsum = sum (FPlistForARead)
		rec = TPsum / (TPsum + FNsum) if (TPsum + FNsum) != 0 else 0
		prec = TPsum / (TPsum + FPsum) if (TPsum + FPsum) != 0 else 0
		if missingInRead != 0:
			missingSize.append(missingInRead)
		corBRate = sum(corBasesForARead)/(sum(corBasesForARead) + sum(uncorBasesForARead))
		outPerReadMetrics.write(str(rec) + " recall\n")
		outPerReadMetrics.write(str(prec) + " precision\n")
		outPerReadMetrics.write(str(corBRate) + " correct_rate\n")
		recall += rec
		precision += prec
		corBasesRate += corBRate
	nbReadsToDivide += 1
	GCRateRef.append(GCRateRefRead)
	GCRateCorr.append(GCRateCorrRead)
	return recall, precision, corBasesRate, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, nbReadsToDivide

def computeMetrics(fileName, outMSAProfile, outPerReadMetrics, correctedFileName):
	#global metrics to return
	precision = 0
	recall = 0
	corBasesRate = 0
	missingSize = []
	smallReadNumber = 0
	GCRateRef = []
	GCRateCorr = []
	indelsubsUncorr = 0
	indelsubsCorr = 0
	#######################
	msa = open(fileName, 'r')
	nbLines = 0
	lines = msa.readlines()
	correctedReadsList = getCorrectedReads(correctedFileName)
	upperCasePositions = getUpperCasePositions(correctedReadsList, lines)
	indelsubsCorr = [0,0,0] #ins del subs
	indelsubsUncorr = [0,0,0]
	headerNo = 0
	readNo = 0
	nbReadsToDivide = 0
	prevHeader = ""
	sameLastHeader = False #useful for the last triplet
	while nbLines < len(lines):
		if not ">" in lines[nbLines]:
			sameLastHeader = False
			reference = lines[nbLines].rstrip() # get msa for ref
			nbLines += 2
			uncorrected =  lines[nbLines].rstrip() # msa for uncorrected
			nbLines += 2
			corrected = lines[nbLines].rstrip() # msa for corrected
			nbLines += 1 #go to next header
			lenCorrected = getLen(corrected)
			lenReference = getLen(reference)
			stretches = findGapStretches(corrected, reference)
			correctedPositions, existingCorrectedReadPositions = getCorrectedPositions(stretches, len(corrected), readNo, upperCasePositions, reference)
			# if the corrected read (or this part of a split read) is long enough
			FP, TP, FN, corBases, uncorBases, GCRateRefRead, GCRateCorrRead, insU, deleU, subsU, insC, deleC, subsC = getTPFNFP(reference, uncorrected, corrected, correctedPositions, existingCorrectedReadPositions)
			# add insertion/deletions/substitutions that have been counted
			indelsubsCorr[0] += insC
			indelsubsCorr[1] += deleC
			indelsubsCorr[2] += subsC
			indelsubsUncorr[0] += insU
			indelsubsUncorr[1] += deleU
			indelsubsUncorr[2] += subsU
			if headerNo == prevHeader: #read in several parts (split) : do not output, only store information
				sameLastHeader = True
				# we add to list the several rates measured on each part
				corBasesForARead.append(corBases)
				uncorBasesForARead.append(uncorBases)
				FPlistForARead.append(FP)
				TPlistForARead.append(TP)
				FNlistForARead.append(FN)
				lenPrevReference = lenReference
				existingCorrectedPositionsInThisRead = [any(tup) for tup in zip(existingCorrectedReadPositions, existingCorrectedPositionsInThisRead)] #logical OR, positions that exist in the corrected read 
				correctedPositionsRead = [any(tup) for tup in zip(correctedPositionsRead, correctedPositions)] #logical OR, positions that are corrected in the read
			else: # end of previous read in several parts or end of previous simple read or first triplet => output for previous read and start to store info for current read
				if prevHeader != "":
					# output info for previous read
					recall, precision, corBasesRate, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, nbReadsToDivide = outputMetrics(recall, precision, corBasesRate, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, lenPrevReference, existingCorrectedPositionsInThisRead, FPlistForARead, TPlistForARead, FNlistForARead, corBasesForARead,uncorBasesForARead, GCRateRefRead, GCRateCorrRead, nbReadsToDivide)
				# store new info
				corBasesForARead = [corBases]
				uncorBasesForARead = [uncorBases]
				FPlistForARead = [FP]
				TPlistForARead = [TP]
				FNlistForARead = [FN]
				lenPrevReference = lenReference
				existingCorrectedPositionsInThisRead = existingCorrectedReadPositions
				correctedPositionsRead = correctedPositions
				prevHeader = headerNo
				if len(lines) == 6:
					recall, precision, corBasesRate, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, nbReadsToDivide = outputMetrics(recall, precision, corBasesRate, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, lenPrevReference, existingCorrectedPositionsInThisRead, FPlistForARead, TPlistForARead, FNlistForARead, corBasesForARead,  uncorBasesForARead, GCRateRefRead, GCRateCorrRead, nbReadsToDivide)
			readNo += 1
				
		else:
			headerNo = lines[nbLines].split(">")[1].split(" ")[0]
			nbLines += 1
	if sameLastHeader: # we must output info for the last read
		recall, precision, corBasesRate, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, nbReadsToDivide = outputMetrics(recall, precision, corBasesRate, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, lenPrevReference, existingCorrectedPositionsInThisRead, FPlistForARead, TPlistForARead, FNlistForARead, corBasesForARead, uncorBasesForARead, GCRateRefRead, GCRateCorrRead, nbReadsToDivide)
	# compute global metrics
	GCRateRef = sum(GCRateRef) / len(GCRateRef)
	GCRateCorr = sum(GCRateCorr) / len(GCRateCorr)
	recall = recall*1.0 / nbReadsToDivide if nbReadsToDivide != 0 else 0
	precision = precision*1.0 / nbReadsToDivide if nbReadsToDivide != 0 else 0
	corBasesRate = corBasesRate*1.0 / nbReadsToDivide if nbReadsToDivide != 0 else 0
	return precision, recall, corBasesRate, missingSize,  GCRateRef, GCRateCorr, indelsubsUncorr, indelsubsCorr


# get the position of nt in uppercase to compute recall and precision only at these positions
def getUpperCasePositions(correctedReadsList, lines):
	upperCasePositions = [] # positions to take into account in the msa
	nbLines = 5 # starting at the first corrected sequence
	headerNo = lines[0].split(">")[1].split(" ")[0]
	headerPrev = ""
	while nbLines < len(lines):
		headerNo = lines[nbLines - 1].split(">")[1].split(" ")[0]
		correctedMsa = lines[nbLines].rstrip()
		#TOVERIFY
		if headerNo == headerPrev:
			index += 1
		else:
			index = 0
		if headerNo in correctedReadsList.keys():
			correctedReadSequence = correctedReadsList[headerNo][index]
			headerPrev = headerNo
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
		else:
			upperCasePositions[-1] = [False] * len(correctedMsa)
	return upperCasePositions


# add to the uppercase positions the positions where there is no stretch of "." , i.e. all positions where recall and precision are actually computed
def getCorrectedPositions(stretches, msaLineLen, readNo, upperCasePositions, reference):
	correctedPositions = copy.copy(upperCasePositions[readNo]) #all positions in upper case in the corrected read
	existingCorrectedPositions = [True] * len(correctedPositions)
	positionsToRemove = list()
	if len(stretches.keys()) > 0:  # split read (or trimmed)
		for pos in stretches.keys(): 
			positionsToRemove.append([pos, stretches[pos]]) #interval(s)) in which the corrected read sequence does not exist
		for interv in positionsToRemove:
			for i in range(interv[0], interv[1] ):
				existingCorrectedPositions[i] = False # remove regions where there is no corrected sequence
				correctedPositions[i] = False # remove regions where there is no corrected sequence (split/trimmed) from corrected regions
	return correctedPositions, existingCorrectedPositions
