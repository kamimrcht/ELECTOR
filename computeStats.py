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
import statistics

#SIZE_CORRECTED_READ_THRESHOLD = 0.1

#THRESH = 100
#THRESH2 = 400

THRESH=5
THRESH2=20

#@Camille, modif ici pour les LR corrigés étendus, avec d'autres fonctions de détection de gaps
def nbLeftGaps(sequence):
	nbGaps = 0
	nbNt = 0
	totalGaps = 0
	i = 0
	while (i < len(sequence) and nbNt <= THRESH):
		if sequence[i] == ".":
			nbGaps += 1
			nbNt = 0
		else:
			if (nbGaps >= 5):
				totalGaps = i
			nbGaps = 0
			nbNt += 1
		i += 1

	return totalGaps

#@Camille, modif ici pour les LR corrigés étendus, avec d'autres fonctions de détection de gaps
def nbRightGaps(sequence):
	nbGaps = 0
	nbNt = 0
	totalGaps = 0
	i = len(sequence) - 1
	while (i >= 0 and nbNt <= THRESH):
		if sequence[i] == ".":
			nbGaps += 1
			nbNt = 0
		else:
			if (nbGaps >= 5):
				totalGaps = len(sequence) - i
			nbGaps = 0
			nbNt += 1
		i -= 1

	return totalGaps

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
#@Camille j'ai aussi fait un peu de ménage ici
def outputRecallPrecision( correctedFileName, outDir, logFile, smallReadNumber, wronglyCorrectedReadsNumber, reportedHomopolThreshold, SIZE_CORRECTED_READ_THRESHOLD,  fileSizeName,beg=0, end=0, soft=None):
	if soft is not None:
		#recall, precision
		outMetrics = open(outDir + "/" + soft + "_per_read_metrics.txt", 'w')
		outMetrics.write("score metric\n")
		nbReads, throughput, precision, recall, corBasesRate, errorRate, extendedBases, missingSize,  GCRateRef, GCRateCorr,  indelsubsUncorr, indelsubsCorr,  trimmedOrSplit, ratioHomopolymers, lenAllCorrectedReads, globalRec, globalPrec  = computeMetrics(outDir + "/msa_" + soft + ".fa", outMetrics, correctedFileName, reportedHomopolThreshold )

	else:
		outMetrics = open(outDir + "/per_read_metrics.txt", 'w')
		outMetrics.write("score metric\n")
		nbReads, throughput, precision, recall, corBasesRate, errorRate, extendedBases, missingSize,  GCRateRef, GCRateCorr, indelsubsUncorr, indelsubsCorr, trimmedOrSplit, ratioHomopolymers, lenAllCorrectedReads, globalRec, globalPrec  = computeMetrics(outDir + "/msa.fa", outMetrics, correctedFileName, reportedHomopolThreshold)

	# read lengths
	outputReadSizeDistribution(correctedFileName, fileSizeName, outDir, trimmedOrSplit, lenAllCorrectedReads)

	outMetrics.close()
	meanMissingSize = 0
	if len(missingSize) > 0:
		meanMissingSize = round(sum(missingSize)/len(missingSize),1)
	meanExtendedBases = 0
	if len(extendedBases) > 0:
		meanExtendedBases = round(sum(extendedBases)/len(extendedBases),1)
	if soft is not None:
		print(soft)

	

	recall = round(recall, 7)
	precision = round(precision, 7)
	corBasesRate = round(corBasesRate, 7)
	errorRate = round(errorRate, 7)
	GCRateRef = round(GCRateRef * 100, 7)
	GCRateCorr = round(GCRateCorr * 100, 7)

	globalRecall =  round(globalRec, 7)
	globalPrecision =  round(globalPrec, 7)
	print("*********** SUMMARY ***********")
	print("Assessed reads: ", str(nbReads))
	print("Throughput: ", str(throughput))
	print("Recall (computed only on corrected bases):", str(recall))
	print("Global recall (computed on all bases):", str(globalRecall))
	print("Precision (computed only on corrected bases):", str(precision))
	print("Global precision (computed on all bases:", str(globalPrecision))
	print("Average correct bases rate:", str(corBasesRate))
	print("Overall error rate: ", str(errorRate))
	print("Number of trimmed/split reads:" , str(trimmedOrSplit))
	print("Mean missing size in trimmed/split reads:" , str(meanMissingSize))
	print("Number of over-corrected reads by extention: ", str(len(extendedBases)))
	print("Mean extension size in over-corrected reads: ", str(meanExtendedBases))
	print("%GC in reference reads: ", str(GCRateRef))
	print("%GC in corrected reads: ", str(GCRateCorr))
	print("Number of corrected reads which length is <", SIZE_CORRECTED_READ_THRESHOLD*100,"% of the original read:", smallReadNumber)
	print("Number of very low quality corrected reads: ", wronglyCorrectedReadsNumber)
	print("Number of insertions in uncorrected: ", str(indelsubsUncorr[0]))
	print("Number of insertions in corrected: ", str(indelsubsCorr[0]))
	print("Number of deletions in uncorrected: ", str(indelsubsUncorr[1]))
	print("Number of deletions in corrected: ", str(indelsubsCorr[1]))
	print("Number of substitutions in uncorrected: ", str(indelsubsUncorr[2]))
	print("Number of substitutions in corrected: ", str(indelsubsCorr[2]))
	#TODO
	print("Ratio of homopolymer sizes in corrected vs reference:", str(ratioHomopolymers))
	
	logFile.write("*********** SUMMARY ***********\n" + "Assessed reads: " + str(nbReads) +"\nThroughput: " + str(throughput) + "\nRecall (computed only on corrected bases):" + str(recall) + "\nGlobal recall (computed on all bases):" + str(globalRecall) + "\nPrecision (computed only on corrected bases):" + str(precision) + "\nGlobal precision (computed on all bases:" + str(globalPrecision) + "\nAverage correct bases rate:" + str(round(corBasesRate,5)) + "\nOverall error rate: " + str(errorRate) + "\nNumber of trimmed/split reads:" + str(trimmedOrSplit) + "\nMean missing size in trimmed/split reads:" + str(meanMissingSize) + "\nNumber of over-corrected reads by extention: " + str(len(extendedBases)) + "\nMean extension size in over-corrected reads: " + str(meanExtendedBases) + "\n%GC in reference reads: " + str(GCRateRef * 100) + "\n%GC in corrected reads: " + str(GCRateCorr * 100) + "\nNumber of corrected reads which length is <" + str(SIZE_CORRECTED_READ_THRESHOLD*100) + "% of the original read:" + str(smallReadNumber) + "\nNumber of very low quality corrected reads: " + str(wronglyCorrectedReadsNumber) + "\nNumber of insertions in uncorrected: " + str(indelsubsUncorr[0]) +"\nNumber of insertions in corrected: " + str(indelsubsCorr[0]) + "\nNumber of deletions in uncorrected: " + str(indelsubsUncorr[1]) +"\nNumber of deletions in corrected: " + str(indelsubsCorr[1]) + "\nNumber of substitutions in uncorrected: " + str(indelsubsUncorr[2]) +"\nNumber of substitutions in corrected: " + str(indelsubsCorr[2]) + 	"Ratio of homopolymer sizes in corrected vs reference: " + str(ratioHomopolymers) + "\n")
	return nbReads, throughput, precision, recall, corBasesRate, errorRate, smallReadNumber, wronglyCorrectedReadsNumber, GCRateRef, GCRateCorr, str(len(missingSize)) , meanMissingSize, str(len(extendedBases)), meanExtendedBases, SIZE_CORRECTED_READ_THRESHOLD,  indelsubsUncorr, indelsubsCorr, trimmedOrSplit, ratioHomopolymers, globalRecall, globalPrecision


def getLen(sequenceMsa):
	return len(sequenceMsa) - sequenceMsa.count('.')

# Compute the length distribution of uncorrected and corrected reads
def outputReadSizeDistribution(correctedFileName, outFileName, outDir, trimmedOrSplit, lenAllReads):
	out = open(outDir + "/" + outFileName, 'w')
	out.write("size type\n")
	for readSize in lenAllReads:
		out.write(str(readSize) + " reads\n")
	if trimmedOrSplit != 0:
	#~ unco = open(uncorrectedFileName)
		cor = open(correctedFileName)
		
		
		#~ l = unco.readline()
		#~ while l != "":
			#~ l = unco.readline()[:-1]
			#~ out.write(str(len(l)) + " uncorrected\n")
			#~ l = unco.readline()
		l = cor.readline()
		while l != "":
			l = cor.readline()[:-1]
			out.write(str(len(l)) + " sequences\n")
			l = cor.readline()
		#~ unco.close()
		cor.close()
	out.close()


# compute ins, del, subs
def indels(ntRef, ntUnco, ntResult,  existingCorrectedPositions, position, insU, deleU, subsU, insC, deleC, subsC, reported, detectedHomopolymer, endOfHomopolRef, okToReportRef, reportedThreshold):
	#compute indels in uncorrected reads
	# and homopolymers
	endOfHomopolResult = True
	okToAppendR = False
	okToAppendC = False
	##### homopolymers in ref ######
	if existingCorrectedPositions[position]:
		if ntRef != '.':
			if ntRef == reported[0][-1]:
				okToAppendR = True #elongation of the repeated chain
				if len(reported[0]) + 1 >= reportedThreshold: #repeated chain of length > THRESHOLD = a homopolymer
					okToReportRef = True #we can report the homopolymer from now, however it can still be elongate
			else:
				if okToReportRef: #we were in a homopolymer and it just terminated at this base => report it
					endOfHomopolRef = True
	###############################
		if ntUnco != ntRef:
			if ntRef == ".":
				insU += 1
			if ntUnco != "." :
				subsU += 1
			else:
				deleU += 1
	#compute only indels in parts of the MSA that actually correspond to a portion that exist in the corrected read
		if ntResult != ntRef:
			if ntRef == ".":
				insC += 1
			else:
				if ntResult != "." :
					subsC += 1
				else:
					deleC += 1

		##########" homopolymer in corrected ############
		if ntResult != '.':
			if ntResult == reported[1][-1]:
				okToAppendC = True #elongation of the repeated chain
				endOfHomopolResult = False #we still are elongating of the repeated chain for corrected read => do not stop
			#~ else:
				#~ endOfHomopolResult = False
	if okToAppendC or okToAppendR:
		reported[0].append(ntRef)
		reported[1].append(ntResult)
	else:
		if not (endOfHomopolRef and endOfHomopolResult):
			if (not endOfHomopolRef and ntRef != "."):
				reported= [[ntRef],[ntResult]]
		
	if endOfHomopolRef and endOfHomopolResult:  #finally report an homopolymer
		detectedHomopolymer = True
		ntH = max(set(reported[0]), key=reported[0].count)
		if ntH == '.':
			filtered = list(filter(lambda x:  x!= '.', reported[0]))
			ntH = max(set(filtered), key=filtered.count)
		sizeHomopolR = [0]
		sizeHomopolC = [0]
		for nt, nt2 in zip(reported[0], reported[1]):
			if nt == ntH:
				sizeHomopolR[-1] += 1
			elif nt != ".":
				sizeHomopolR.append(0)
			if nt2 == ntH:
				sizeHomopolC[-1] += 1
			elif nt2 != ".":
				sizeHomopolC.append(0)
		sizeHomopolReference = max(sizeHomopolR)
		sizeHomopolCorrected = max(sizeHomopolC)
		reported = [sizeHomopolReference, sizeHomopolCorrected]
				
		
	return insU, deleU, subsU, insC, deleC, subsC, reported, detectedHomopolymer, endOfHomopolRef, okToReportRef

# compute fp, tp, fn
def getCorrectionAtEachPosition(ntRef, ntUnco, ntResult, correctedPositions, existingCorrectedPositions, position, corBases, uncorBase, FP, FN, TP, globalFP, globalFN, globalTP):
	#compute FN,FP,TP only on corrected parts
	if correctedPositions[position]:
		if ntRef == ntUnco:  #no error
			if ntUnco != ntResult: #FP
				FP += 1
				globalFP += 1
				uncorBase += 1
				corBases -= 1
			# else good nt not corrected = ok
			#else:
			#	corBases +=1
		else: #error
			if ntRef == ntResult: #error corrected
				TP += 1
				globalTP += 1
			#	corBases += 1
			else:
				if ntUnco == ntResult: # error not corrected
					FN += 1
					globalFN += 1
				else: #new error introduced by corrector
					FP += 1
					globalFP += 1
				uncorBase += 1
				corBases -= 1
	#correct base rate is computed everywhere
	else:
		if existingCorrectedPositions[position]:
			#if ntRef == ntResult:
			#	corBases += 1
			#else:
			if ntRef != ntResult:
				uncorBase += 1
				corBases -= 1
			if ntRef == ntUnco:  #no error
				if ntUnco != ntResult: #FP
					globalFP += 1
			else: #error
				if ntRef == ntResult: #error corrected
					globalTP += 1
				#	corBases += 1
				else:
					if ntUnco == ntResult: # error not corrected
						globalFN += 1
					else: #new error introduced by corrector
						globalFP += 1
	return corBases, uncorBase, FP, FN, TP, globalFP, globalFN, globalTP


# get insertion deletion substitution FP, FN, TP and GC rates for a triplet
def getTPFNFP(reference, corrected, uncorrected,  correctedPositions, existingCorrectedPositions, reportedThreshold, ratioHomopolymers):
	position = 0
	corBases = getLen(corrected)
	uncorBases = 0
	FN = 0
	FP = 0
	TP = 0
	globalFN = 0
	globalFP = 0
	globalTP = 0
	GCSumRef = 0
	GCSumCorr = 0
	insU = 0
	deleU = 0
	subsU = 0
	insC = 0
	deleC = 0
	subsC = 0
	detectedHomopolymer = False
	okToReportRef = False
	endOfHomopolRef = False
	reported = [['x'],['x']]
	for ntRef, ntResult, ntUnco in zip(reference, corrected, uncorrected): #zip(reference, uncorrected, corrected):
		if ntRef.upper() == "G" or ntRef.upper() == "C":
			GCSumRef += 1
		if ntResult.upper() ==  "G" or ntResult.upper() == "C":
			GCSumCorr += 1
		#insertion deletion substitution
		insU, deleU, subsU, insC, deleC, subsC, reported, detectedHomopolymer, endOfHomopolRef, okToReportRef= indels(ntRef, ntUnco, ntResult,  existingCorrectedPositions, position, insU, deleU, subsU, insC, deleC, subsC, reported, detectedHomopolymer, endOfHomopolRef, okToReportRef, reportedThreshold)

		if detectedHomopolymer:
			detectedHomopolymer = False
			okToReportRef = False
			endOfHomopolRef = False
			ratioHomopolymers.append(round(reported[1]*1.0/reported[0],2))
			reported = [ [ntRef], [ntResult]]
		# HERE replace homopol in U by those in R and if homopol is reported in R dividec hC/hR and append a vector (this vector can be the same for all reads)
		#FP, FN, TP
		corBases, uncorBases, FP, FN, TP,globalFP, globalFN, globalTP = getCorrectionAtEachPosition(ntRef, ntUnco, ntResult, correctedPositions,  existingCorrectedPositions, position,  corBases, uncorBases, FP, FN, TP,globalFP, globalFN, globalTP)
		#print(str(corBases))
		position += 1
	GCRateRef = round(GCSumRef * 1.0 / getLen(reference),3)
	GCRateCorr = round(GCSumCorr * 1.0 / getLen(corrected),3)
	#~ return FP, TP, FN, corBases, uncorBases, GCRateRef, GCRateCorr, insU, deleU, subsU, insC, deleC, subsC, homopolymersUInser, homopolymersUDele, homopolymersCInser, homopolymersCDele
	return FP, TP, FN, corBases, uncorBases, GCRateRef, GCRateCorr, insU, deleU, subsU, insC, deleC, subsC, ratioHomopolymers, globalFP, globalFN, globalTP


def outputMetrics(recall, precision, corBasesRate, missingSize, extendedBases, totalGaps, GCRateRef, GCRateCorr, outPerReadMetrics, lenReference, existingCorrectedPositionsInThisRead, FPlistForARead, TPlistForARead, FNlistForARead, corBasesForARead, uncorBasesForARead, GCRateRefRead, GCRateCorrRead, nbReadsToDivide, totalCorBases, totalUncorBases, sameLastHeader,lenCorrected, lenAllCorrectedReads,globalFPlistForARead, globalTPlistForARead, globalFNlistForARead, trimmedOrSplit):
	missingInRead = existingCorrectedPositionsInThisRead.count(False) # final sum
	if not sameLastHeader and missingInRead != 0:
		trimmedOrSplit += 1
	#~ else:
		#~ missingInRead = 0
	if FPlistForARead != [] or TPlistForARead != [] or FNlistForARead != [] or corBasesForARead != []:
		FNsum = sum(FNlistForARead)
		TPsum = sum(TPlistForARead)
		FPsum = sum (FPlistForARead)
		rec = TPsum / (TPsum + FNsum) if (TPsum + FNsum) != 0 else 0
		prec = TPsum / (TPsum + FPsum) if (TPsum + FPsum) != 0 else 0
		if missingInRead != 0:
			missingSize.append(missingInRead)
		if totalGaps > 0:
			extendedBases.append(totalGaps)
		corBRate = sum(corBasesForARead)/(sum(corBasesForARead) + sum(uncorBasesForARead))
		outPerReadMetrics.write(str(rec) + " recall\n")
		outPerReadMetrics.write(str(prec) + " precision\n")
		outPerReadMetrics.write(str(corBRate) + " correct_rate\n")
		recall.append(rec)
		precision.append(prec)
		corBasesRate.append(corBRate)
		totalCorBases += sum(corBasesForARead)
		totalUncorBases += sum(uncorBasesForARead)
	nbReadsToDivide += 1
	if globalFPlistForARead != [] or globalTPlistForARead != [] or globalFNlistForARead != [] :
		globalFNsum = sum(globalFNlistForARead) + missingInRead
		globalTPsum = sum(globalTPlistForARead)
		globalFPsum = sum (globalFPlistForARead)
		globalRec = globalTPsum / (globalTPsum + globalFNsum) if (globalTPsum + globalFNsum) != 0 else 0
		globalPrec = globalTPsum / (globalTPsum + globalFPsum) if (globalTPsum + globalFPsum) != 0 else 0
	GCRateRef.append(GCRateRefRead)
	GCRateCorr.append(GCRateCorrRead)
	lenAllCorrectedReads.append(lenCorrected)
	return recall, precision, corBasesRate, missingSize, extendedBases, GCRateRef, GCRateCorr, outPerReadMetrics, nbReadsToDivide, totalCorBases, totalUncorBases, lenAllCorrectedReads, globalRec, globalPrec, trimmedOrSplit

def computeMetrics(fileName, outPerReadMetrics, correctedFileName, reportedThreshold):
	#global metrics to return
	precision = []
	recall = []
	corBasesRate = []
	totalCorBases = 0
	totalUncorBases = 0
	missingSize = []
	extendedBases = []
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
	numberHomopolymersInserInCorrected = 0
	numberHomopolymersDeleInCorrected = 0
	numberHomopolymersInserInUncorrected = 0
	numberHomopolymersDeleInUncorrected = 0
	meanLengthDeleHomopolymersInUncorrected = []
	meanLengthInserHomopolymersInUncorrected = []
	meanLengthInserHomopolymersInCorrected = []
	meanLengthDeleHomopolymersInCorrected = []
	trimmedOrSplit = 0
	ratioHomopolymers = []
	lenCorrected = 0
	lenAllCorrectedReads = [] #this length is computing by adding the sizes of split reads if needed
	while nbLines < len(lines):
		if not ">" in lines[nbLines]:
			if sameLastHeader:
				trimmedOrSplit += 1
			#~ sameLastHeader = False
			reference = lines[nbLines].rstrip() # get msa for ref
			nbLines += 2
			corrected =  lines[nbLines].rstrip() # msa for uncorrected
			nbLines += 2
			uncorrected = lines[nbLines].rstrip() # msa for corrected

			#@Camille, modif ici pour les LR corrigés étendus, avec d'autres fonctions de détection de gaps
			#Strip "." at the beginning and at the end of reference / uncorrected => caused by an extended corrected read
			#Also strip nts from the corrected sequence, to allow computation of the next step, and record them to output them
			refGapsLeft = nbLeftGaps(reference)
			uncoGapsLeft = nbLeftGaps(uncorrected)
			gapsLeft = min(refGapsLeft, uncoGapsLeft);
			if (gapsLeft >= THRESH2):
				reference = reference[gapsLeft:-1]
				uncorrected = uncorrected[gapsLeft:-1]
				corrected = corrected[gapsLeft:-1]

			refGapsRight = nbRightGaps(reference)
			uncoGapsRight = nbRightGaps(uncorrected)
			gapsRight = min(refGapsRight, uncoGapsRight)
			if (gapsRight >= THRESH2):
				reference = reference[0:len(reference) - gapsRight]
				uncorrected = uncorrected[0:len(uncorrected) - gapsRight]
				corrected = corrected[0:len(corrected) - gapsRight]
			
			nbLines += 1 #go to next header
			#~ lenCorrected = getLen(corrected)
			lenReference = getLen(reference)
			stretches = findGapStretches(corrected, reference)
			correctedPositionsRead, existingCorrectedPositionsInThisRead = getCorrectedPositions(stretches, len(corrected), readNo, upperCasePositions, reference)
			
			FP, TP, FN, corBases, uncorBases, GCRateRefRead, GCRateCorrRead, insU, deleU, subsU, insC, deleC, subsC, ratioHomopolymers,globalFP, globalFN, globalTP = getTPFNFP(reference, corrected, uncorrected, correctedPositionsRead, existingCorrectedPositionsInThisRead, reportedThreshold, ratioHomopolymers)

			# add insertion/deletions/substitutions that have been counted
			indelsubsCorr[0] += insC
			indelsubsCorr[1] += deleC
			indelsubsCorr[2] += subsC
			indelsubsUncorr[0] += insU
			indelsubsUncorr[1] += deleU
			indelsubsUncorr[2] += subsU
			# print(headerNo)
			# print(prevHeader)
			if headerNo == prevHeader: #read in several parts (split) : do not output, only store information
				sameLastHeader = True
				# we add to list the several rates measured on each part
				#~ print("haha this is why")
				corBasesForARead.append(corBases)
				uncorBasesForARead.append(uncorBases)
				FPlistForARead.append(FP)
				TPlistForARead.append(TP)
				FNlistForARead.append(FN)
				globalFPlistForARead.append(globalFP)
				globalTPlistForARead.append(globalTP)
				globalFNlistForARead.append(globalFN)
				totalGaps = totalGaps + gapsLeft + gapsRight
				lenPrevReference = lenReference
				existingCorrectedPositionsInThisRead = [any(tup) for tup in zip(existingCorrectedReadPositions, existingCorrectedPositionsInThisRead)] #logical OR, positions that exist in the corrected read 
				correctedPositionsRead = [any(tup) for tup in zip(correctedPositionsRead, correctedPositions)] #logical OR, positions that are corrected in the read
				lenRead += getLen(corrected)
			else: # end of previous read in several parts or end of previous simple read or first triplet => output for previous read and start to store info for current read
					
				if prevHeader != "":
					# output info for previous read
					recall, precision, corBasesRate, missingSize, extendedBases, GCRateRef, GCRateCorr, outPerReadMetrics, nbReadsToDivide, totalCorBases, totalUncorBases, lenAllCorrectedReads, globalRec, globalPrec, trimmedOrSplit = outputMetrics(recall, precision, corBasesRate, missingSize, extendedBases, totalGaps, GCRateRef, GCRateCorr, outPerReadMetrics, lenPrevReference, existingCorrectedPositionsInThisRead, FPlistForARead, TPlistForARead, FNlistForARead, corBasesForARead,uncorBasesForARead, GCRateRefRead, GCRateCorrRead, nbReadsToDivide, totalCorBases, totalUncorBases, sameLastHeader, lenCorrected, lenAllCorrectedReads,globalFPlistForARead, globalTPlistForARead, globalFNlistForARead, trimmedOrSplit)
				# store new info
				lenCorrected = getLen(corrected)
				sameLastHeader = False
				corBasesForARead = [corBases]
				uncorBasesForARead = [uncorBases]
				FPlistForARead = [FP]
				TPlistForARead = [TP]
				FNlistForARead = [FN]
				globalFPlistForARead = [globalFP]
				globalTPlistForARead = [globalTP]
				globalFNlistForARead = [globalFN]
				totalGaps = gapsLeft + gapsRight
				lenPrevReference = lenReference
				existingCorrectedReadPositions = existingCorrectedPositionsInThisRead
				correctedPositions = correctedPositionsRead
				prevHeader = headerNo
				if len(lines) == 6:
					recall, precision, corBasesRate, missingSize, extendedBases, GCRateRef, GCRateCorr, outPerReadMetrics, nbReadsToDivide, totalCorBases, totalUncorBases, lenAllCorrectedReads, globalRec, globalPrec, trimmedOrSplit = outputMetrics(recall, precision, corBasesRate, missingSize, extendedBases, totalGaps, GCRateRef, GCRateCorr, outPerReadMetrics, lenPrevReference, existingCorrectedPositionsInThisRead, FPlistForARead, TPlistForARead, FNlistForARead, corBasesForARead,  uncorBasesForARead, GCRateRefRead, GCRateCorrRead, nbReadsToDivide, totalCorBases, totalUncorBases, sameLastHeader, lenCorrected, lenAllCorrectedReads, globalFPlistForARead, globalTPlistForARead, globalFNlistForARead, trimmedOrSplit)
			readNo += 1
				
		else:
			headerNo = lines[nbLines].split(">")[1].split(" ")[0]
			nbLines += 1
	if sameLastHeader or prevHeader != "": # we must output info for the last read
		if sameLastHeader:
			trimmedOrSplit += 1
			lenCorrected += getLen(corrected)
		lenCorrected = getLen(corrected)
		recall, precision, corBasesRate, missingSize, extendedBases, GCRateRef, GCRateCorr, outPerReadMetrics, nbReadsToDivide, totalCorBases, totalUncorBases, lenAllCorrectedReads, globalRec, globalPrec, trimmedOrSplit = outputMetrics(recall, precision, corBasesRate, missingSize, extendedBases, totalGaps, GCRateRef, GCRateCorr, outPerReadMetrics, lenPrevReference, existingCorrectedPositionsInThisRead, FPlistForARead, TPlistForARead, FNlistForARead, corBasesForARead, uncorBasesForARead, GCRateRefRead, GCRateCorrRead, nbReadsToDivide, totalCorBases, totalUncorBases, sameLastHeader, lenCorrected, lenAllCorrectedReads,globalFPlistForARead, globalTPlistForARead, globalFNlistForARead, trimmedOrSplit)
	# compute global metrics
	GCRateRef = round(sum(GCRateRef) / len(GCRateRef),3)
	GCRateCorr = round(sum(GCRateCorr) / len(GCRateCorr),3)
	recall = sum(recall)*1.0 / nbReadsToDivide if nbReadsToDivide != 0 else 0
	precision = sum(precision)*1.0 / nbReadsToDivide if nbReadsToDivide != 0 else 0
	corBasesRate = sum(corBasesRate)*1.0 / nbReadsToDivide if nbReadsToDivide != 0 else 0
	throughput = totalCorBases + totalUncorBases
	errorRate = 1 - (totalCorBases / throughput)
	if len(ratioHomopolymers) > 1:
		meanRatioHomopolymers = statistics.mean(ratioHomopolymers)
	else:
		meanRatioHomopolymers = 1
	return nbReadsToDivide, throughput, precision, recall, corBasesRate, errorRate, extendedBases, missingSize,  GCRateRef, GCRateCorr, indelsubsUncorr, indelsubsCorr, trimmedOrSplit, meanRatioHomopolymers, lenAllCorrectedReads, globalRec, globalPrec


# get the position of nt in uppercase to compute recall and precision only at these positions
def getUpperCasePositions(correctedReadsList, lines):
	upperCasePositions = [] # positions to take into account in the msa
	nbLines = 3 # starting at the first corrected sequence
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
