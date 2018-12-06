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
import utils
#~ import ast
#SIZE_CORRECTED_READ_THRESHOLD = 0.1

#THRESH = 100
#THRESH2 = 400

THRESH=5
THRESH2=20



def getSplit(fileName):
	cmdGrep = """grep ">" """ + fileName + """ | uniq -c """
	lines = subprocess.check_output(['bash','-c', cmdGrep]).decode('ascii')
	lines = lines.split('\n')
	#~ lines = subprocess.check_output(['bash','-c', cmdGrep]).decode('ascii').split['\n']
	readToSplit = dict()
	for l in lines:
		if len(l) > 0:
			nb = int(l.split('>')[0].replace(' ','').replace('\t',''))
			read = l.split('>')[1].replace(' ','').replace('\t','')
			readToSplit[read] = int(nb/3)
	return readToSplit

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
			if (nbGaps >= THRESH):
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
			if (nbGaps >= THRESH):
				totalGaps = len(sequence) - i
			nbGaps = 0
			nbNt += 1
		i -= 1

	return totalGaps

# store corrected reads in a list of tuples [(header, sequence)]
#~ def getCorrectedReads(correctedReadsFileName):
	#~ fastaTuple = dict()
	#~ with open(correctedReadsFileName) as fileIn:
		#~ for line in fileIn:
			#~ if ">" in line:
				#~ header = line.rstrip().split(' ')[0][1:]
				
			#~ else:
				#~ read = line.rstrip()
				#~ if read != "":
					
					#~ if header not in fastaTuple.keys():
						#~ fastaTuple[header] = [read]
					#~ else:
						#~ fastaTuple[header].append(read)
	#~ return fastaTuple


# find long stretches of "." = trimmed or split reads, and return the coordinates of these regions for the corrected line of a msa
# if reference is not a "." too, else it means it is a gap opened by the uncorrected part
def findGapStretches(correctedSequence, referenceSequence, gapsPositions):
	prev = None
	countGap = 0
	countGapRef = 0
	positionsStretch = []
	pos = 0

	#~ for pos,ntResult in enumerate(correctedSequence):   # look for gaps in splitted/trimmed corrected read
	for ntRef,ntResult in zip(referenceSequence, correctedSequence):   # look for gaps in splitted/trimmed corrected read
		#~ if pos not in gapsPositions:
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
				if countGapRef < THRESH2: #if countGapRef>=THRESH2 it means that a gap is openened both in ref and corrected because of the uncorrected seq in the msa, so this is not a trimmed zone
					if len(positionsStretch) == 0:
						positionsStretch.append([pos-THRESH + 1, pos]) # start new stretch of gap with leftmost position
					else:
						if len(positionsStretch[-1]) == 0:
							positionsStretch[-1].extend((pos-THRESH + 1, pos))
						if len(positionsStretch[-1]) == 2:
							positionsStretch[-1][1] = pos # update position
				#~ else:
					#~ print("hop", pos, countGapRef)
					# forNotExisting.append(pos)
			prev = ntResult
			pos += 1

	tmpStretch = []
	# bords
	for i,s in enumerate(positionsStretch):
		if len(positionsStretch) > 1:
			if len(s) > 0:
				#~ if i == 0:
					if s[0] <= THRESH2:
						tmpStretch.append([0, s[1]])
					#~ else:
						#~ tmpStretch.append([s[0], s[1]])
				#~ elif i == len(positionsStretch) - 1:
					if len(correctedSequence) - s[1] <= THRESH2:
						tmpStretch.append([s[0], len(correctedSequence) - 1])
					else:
						tmpStretch.append([s[0], s[1]])
				#~ else:
					#~ tmpStretch.append([s[0], s[1]])
		elif  len(positionsStretch) == 1:
			if len(s) > 0:
				if s[0] <= THRESH2:
					tmpStretch.append([0, s[1]])
				else:
					tmpStretch.append([s[0], s[1]])
				if len(correctedSequence) - s[1] <= THRESH2:
					tmpStretch[-1][1] = len(correctedSequence) - 1
		
	#merge
	tmpStretch2 = []
	merge = False
	for i,s in enumerate(tmpStretch):
		if i < len(tmpStretch) - 1:
			if tmpStretch[i+1][0] - s[1] <= THRESH:
				tmpStretch2.append([s[0], tmpStretch[i+1][1]])
				merge = True
			else:
				tmpStretch2.append([s[0], s[1]])
				merge = False
	if not merge: #add the last one
		if len(tmpStretch) > 0:
			tmpStretch2.append([tmpStretch[-1][0], tmpStretch[-1][1]])
		
	stretch = dict()
	#~ print("tmpStretch2", tmpStretch2)

	#garder seulement aux bords (=split/trimmed)
	for s in tmpStretch2:
		if s[0] == 0:
			if s[1] - s[0] > THRESH2:
				stretch[s[0]] = s[1]
		elif s[1] == len(correctedSequence) - 1:
			if s[1] - s[0] > THRESH2:
				stretch[s[0]] = s[1]
	#~ if len(tmpStretch2) > 0:
		#~ stretch[tmpStretch2[0][0]] = tmpStretch2[0][1]
		#~ if len(tmpStretch2) > 1:
			#~ stretch[tmpStretch2[-1][0]] = tmpStretch2[-1][1]
	return stretch




# compute recall and precision and writes output files
#@Camille j'ai aussi fait un peu de ménage ici
#main function
def outputRecallPrecision(correctedFileName, outDir, logFile, smallReadNumber, wronglyCorrectedReadsNumber, reportedHomopolThreshold, SIZE_CORRECTED_READ_THRESHOLD,  fileSizeName, clipsNb, beg=0, end=0, soft=None):
	if soft is not None:
		#recall, precision
		outMetrics = open(outDir + "/" + soft + "_per_read_metrics.txt", 'w')
		outMetrics.write("score metric\n")
		readsToSplit = getSplit(outDir + "/msa_" + soft + ".fa")
		nbReads, throughput, uncorThroughput, precision, recall, corBasesRate, errorRate, uncorCorBasesRate, uncorErrorRate, missingSize,  GCRateRef, GCRateCorr,  indelsubsUncorr, indelsubsCorr,  ratioHomopolymers, lenAllCorrectedReads,  countReadSplit, countReadTrimmed, countReadExtended, extendedBasesCount  = computeMetrics(outDir + "/msa_" + soft + ".fa", outMetrics, correctedFileName, reportedHomopolThreshold, clipsNb, readsToSplit )

	else:
		outMetrics = open(outDir + "/per_read_metrics.txt", 'w')
		outMetrics.write("score metric\n")
		readsToSplit = getSplit(outDir + "/msa.fa")

		nbReads, throughput, uncorThroughput, precision, recall, corBasesRate, errorRate, uncorCorBasesRate, uncorErrorRate, missingSize,  GCRateRef, GCRateCorr,  indelsubsUncorr, indelsubsCorr,  ratioHomopolymers, lenAllCorrectedReads,  countReadSplit, countReadTrimmed, countReadExtended, extendedBasesCount  = computeMetrics(outDir + "/msa.fa", outMetrics, correctedFileName, reportedHomopolThreshold, clipsNb, readsToSplit)

	# read lengths
	outputReadSizeDistribution(correctedFileName, fileSizeName, outDir, countReadSplit+countReadTrimmed, lenAllCorrectedReads)

	outMetrics.close()
	meanMissingSize = 0
	if countReadSplit + countReadTrimmed > 0:
		#~ meanMissingSize = round(sum(missingSize)/len(missingSize),1)
		meanMissingSize = round(sum(missingSize)/(countReadSplit + countReadTrimmed),1)
	meanExtendedBases = 0
	#~ if len(extendedBases) > 0:
		#~ meanExtendedBases = round(sum(extendedBases)/len(extendedBases),1)
	if countReadExtended > 0:
		#~ print('count', extendedBasesCount)
		meanExtendedBases = round(sum(extendedBasesCount)/countReadExtended,1)
	if soft is not None:
		print(soft)

	

	recall = round(recall, 7)
	precision = round(precision, 7)
	corBasesRate = round(corBasesRate, 7)
	errorRate = round(errorRate, 7)
	GCRateRef = round(GCRateRef * 100, 7)
	GCRateCorr = round(GCRateCorr * 100, 7)

	print("*********** SUMMARY ***********")
	print("Assessed reads: ", str(nbReads))
	print("Throughput (uncorrected)", str(uncorThroughput))
	print("Throughput (corrected): ", str(throughput))
	print("Recall (computed only on corrected bases):", str(recall))
	#~ print("Global recall (computed on all bases):", str(globalRecall))
	print("Precision (computed only on corrected bases):", str(precision))
	#~ print("Global precision (computed on all bases:", str(globalPrecision))
	print("Average correct bases rate (uncorrected): ", str(uncorCorBasesRate))
	print("Error rate (uncorrected):", str(1-uncorCorBasesRate ))
	print("Average correct bases rate (corrected): ", str(corBasesRate))
	print("Error rate (corrected):", str(1-corBasesRate))
	#~ print("Number of trimmed/split reads:" , str(trimmedOrSplit))
	print("Number of trimmed/split reads:" , str(countReadSplit + countReadTrimmed))
	print("Mean missing size in trimmed/split reads:" , str(meanMissingSize))
	print("Number of over-corrected reads by extention: ", str(countReadExtended))
	#~ print("Number of over-corrected reads by extention: ", str(len(extendedBases)))
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
	print("Ratio of homopolymer sizes in corrected vs reference:", str(ratioHomopolymers))
	
	logFile.write("*********** SUMMARY ***********\n" + "Assessed reads: " + str(nbReads) +"\nThroughput (uncorrected): " + str(uncorThroughput) + "\nThroughput (corrected): " + str(throughput) + "\nRecall (computed only on corrected bases):" + str(recall)  + "\nPrecision (computed only on corrected bases):" + str(precision) +  "\nAverage correct bases rate (uncorrected):" + str(uncorCorBasesRate) + "\nError rate (uncorrected): " + str(1-uncorCorBasesRate) +  "\nAverage correct bases rate (corrected):" + str(corBasesRate) + "\nError rate (corrected): " + str(1-corBasesRate) + "\nNumber of trimmed/split reads:" + str(countReadSplit + countReadTrimmed) + "\nMean missing size in trimmed/split reads:" + str(meanMissingSize) + "\nNumber of over-corrected reads by extention: " + str(countReadExtended) + "\nMean extension size in over-corrected reads: " + str(meanExtendedBases) + "\n%GC in reference reads: " + str(GCRateRef) + "\n%GC in corrected reads: " + str(GCRateCorr) + "\nNumber of corrected reads which length is <" + str(SIZE_CORRECTED_READ_THRESHOLD*100) + "% of the original read:" + str(smallReadNumber) + "\nNumber of very low quality corrected reads: " + str(wronglyCorrectedReadsNumber) + "\nNumber of insertions in uncorrected: " + str(indelsubsUncorr[0]) +"\nNumber of insertions in corrected: " + str(indelsubsCorr[0]) + "\nNumber of deletions in uncorrected: " + str(indelsubsUncorr[1]) +"\nNumber of deletions in corrected: " + str(indelsubsCorr[1]) + "\nNumber of substitutions in uncorrected: " + str(indelsubsUncorr[2]) +"\nNumber of substitutions in corrected: " + str(indelsubsCorr[2]) + 	"\nRatio of homopolymer sizes in corrected vs reference: " + str(ratioHomopolymers) + "\n")
	return nbReads, throughput, precision, recall, corBasesRate, 1-corBasesRate, smallReadNumber, wronglyCorrectedReadsNumber, GCRateRef, GCRateCorr, str(countReadSplit + countReadTrimmed) , meanMissingSize,  str(countReadExtended), meanExtendedBases, SIZE_CORRECTED_READ_THRESHOLD,  indelsubsUncorr, indelsubsCorr, countReadSplit + countReadTrimmed, ratioHomopolymers

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
def indels(ntRef, ntUnco, ntResult,  existingCorrectedPositions, position, insU, deleU, subsU, insC, deleC, subsC, reported, detectedHomopolymer, endOfHomopolRef, okToReportRef, reportedThreshold, gapsPositions):
	#compute indels in uncorrected reads
	# and homopolymers
	endOfHomopolResult = True
	okToAppendR = False
	okToAppendC = False
	
					#~ print("del",position, ntRef, ntResult, ntUnco)
	if existingCorrectedPositions[position]:
		##### homopolymers in ref ######
		if ntRef != '.':
			if ntRef == reported[0][-1]:
				okToAppendR = True #elongation of the repeated chain
				if len(reported[0]) + 1 >= reportedThreshold: #repeated chain of length > THRESHOLD = a homopolymer
					okToReportRef = True #we can report the homopolymer from now, however it can still be elongate
			else:
				if okToReportRef: #we were in a homopolymer and it just terminated at this base => report it
					endOfHomopolRef = True
	
	#compute only indels in parts of the MSA that actually correspond to a portion that exist in the corrected read
		if ntResult != ntRef:
			if ntRef == ".":
				insC += 1
			else:
				if ntResult != "." :
					subsC += 1
				else:
					deleC += 1
	##### uncorrected
		if position not in gapsPositions:
			if ntUnco != ntRef:
				if ntRef == ".":
					insU += 1
				else:
					if ntUnco != "." :
						subsU += 1
						#~ print("sub",position, ntRef, ntResult, ntUnco)

					else:
						deleU += 1
				
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
# def getCorrectionAtEachPosition(ntRef, ntUnco, ntResult, correctedPositions, existingCorrectedPositions, position, corBases, uncorBase, FP, FN, TP, globalFP, globalFN, globalTP):
def getCorrectionAtEachPosition(ntRef, ntUnco, ntResult, existingCorrectedPositions, position, corBases, uncorBase, uncorCorBases, uncorUncorBases, FP, FN, TP):
	#compute FN,FP,TP only on corrected parts
	# if correctedPositions[position]:
	if existingCorrectedPositions[position]:
		if ntRef == ntUnco:  #no error
			if ntUnco != ntResult: #FP
				FP += 1
				uncorBase += 1
			else:
				TP += 1
				corBases += 1
			uncorCorBases += 1
		else: #error
			if ntRef == ntResult: #error corrected
				TP += 1
				corBases += 1
			else:
				if ntUnco == ntResult: # error not corrected
					FN += 1
					FP += 1
				uncorBase += 1
			uncorUncorBases += 1
	#correct base rate is computed everywhere
	#~ else:
		#~ if existingCorrectedPositions[position]:
			#~ if ntRef == ntUnco:  
				#~ if ntUnco != ntResult: #FP
					#~ globalFP += 1
					#~ uncorBase += 1
				#~ else:
					#~ globalTP += 1
					#~ corBases += 1
			#~ else: #error
				#~ if ntRef == ntResult: #error corrected
					#~ globalTP += 1
					#~ corBases += 1
				#~ else:
					#~ if ntUnco == ntResult: # error not corrected
						#~ globalFN += 1
					#~ else: #new error introduced by corrector
						#~ globalFP += 1
					#~ uncorBase += 1
		#~ else:
			#~ globalFP += 1
			#~ globalFN += 1
					#~ corBases -= 1
				#~ corBases -= 1
	#~ if corBases < 0 :
		#~ corBases = 0
	#~ if uncorBase < 0:
		#~ uncorBase = 0
	return corBases, uncorBase, uncorCorBases, uncorUncorBases, FP, FN, TP


# get insertion deletion substitution FP, FN, TP and GC rates for a triplet
#~ def getTPFNFP(reference, corrected, uncorrected,  correctedPositions, existingCorrectedPositions, reportedThreshold, ratioHomopolymers, gapsPositions):
def getTPFNFP(reference, corrected, uncorrected,   existingCorrectedPositions, reportedThreshold, ratioHomopolymers, gapsPositions):
	position = 0
	corBases = 0
	uncorBases = 0
	uncorCorBases = 0
	uncorUncorBases = 0
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
		insU, deleU, subsU, insC, deleC, subsC, reported, detectedHomopolymer, endOfHomopolRef, okToReportRef= indels(ntRef, ntUnco, ntResult,  existingCorrectedPositions, position, insU, deleU, subsU, insC, deleC, subsC, reported, detectedHomopolymer, endOfHomopolRef, okToReportRef, reportedThreshold, gapsPositions)
		if detectedHomopolymer:
			detectedHomopolymer = False
			okToReportRef = False
			endOfHomopolRef = False
			ratioHomopolymers.append(round(reported[1]*1.0/reported[0],2))
			reported = [ [ntRef], [ntResult]]
		# HERE replace homopol in U by those in R and if homopol is reported in R dividec hC/hR and append a vector (this vector can be the same for all reads)
		#FP, FN, TP
		corBases, uncorBases, uncorCorBases, uncorUncorBases, FP, FN, TP = getCorrectionAtEachPosition(ntRef, ntUnco, ntResult,   existingCorrectedPositions, position,  corBases, uncorBases, uncorCorBases, uncorUncorBases, FP, FN, TP)
		position += 1

	GCRateRef = round(GCSumRef * 1.0 / getLen(reference),3)
	GCRateCorr = round(GCSumCorr * 1.0 / getLen(corrected),3)
	return FP, TP, FN, corBases, uncorBases, uncorCorBases, uncorUncorBases, GCRateRef, GCRateCorr, insU, deleU, subsU, insC, deleC, subsC, ratioHomopolymers


#~ def outputMetrics(recall, precision, globalRecall, globalPrecision, corBasesRate, missingInRead, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, FPlistForARead, TPlistForARead, FNlistForARead, corBasesForARead, uncorBasesForARead, totalCorBases, totalUncorBases, globalFPlistForARead, globalTPlistForARead, globalFNlistForARead, GCRateRefRead, GCRateCorrRead):
def outputMetrics(recall, precision,  corBasesRate, uncorCorBasesRate, missingInRead, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, FPlistForARead, TPlistForARead, FNlistForARead, corBasesForARead, uncorBasesForARead, uncorCorBasesForARead, uncorUncorBasesForARead, totalCorBases, totalUncorBases, uncorTotalCorBases, uncorTOtalUncorBases,  GCRateRefRead, GCRateCorrRead):
	if FPlistForARead != [] or TPlistForARead != [] or FNlistForARead != [] or corBasesForARead != []:
		FNsum = sum(FNlistForARead)
		TPsum = sum(TPlistForARead)
		FPsum = sum (FPlistForARead)
		rec = TPsum / (TPsum + FNsum) if (TPsum + FNsum) != 0 else 0
		prec = TPsum / (TPsum + FPsum) if (TPsum + FPsum) != 0 else 0
		if missingInRead != 0:
			missingSize.append(missingInRead)
		corBRate = sum(corBasesForARead)/(sum(corBasesForARead) + sum(uncorBasesForARead))
		uncorCorBRate = sum(uncorCorBasesForARead)/(sum(uncorCorBasesForARead) + sum(uncorUncorBasesForARead))
		outPerReadMetrics.write(str(rec) + " recall\n")
		outPerReadMetrics.write(str(prec) + " precision\n")
		outPerReadMetrics.write(str(corBRate) + " correct_rate\n")
		recall.append(rec)
		precision.append(prec)
		corBasesRate.append(corBRate)
		uncorCorBasesRate.append(uncorCorBRate)
		totalCorBases += sum(corBasesForARead)
		totalUncorBases += sum(uncorBasesForARead)
		uncorTotalCorBases += sum(uncorCorBasesForARead)
		uncorTOtalUncorBases += sum(uncorUncorBasesForARead)
	#~ if globalFPlistForARead != [] or globalTPlistForARead != [] or globalFNlistForARead != [] :
		#~ globalFNsum = sum(globalFNlistForARead) + missingInRead + FNsum
		#~ globalFNsum = sum(globalFNlistForARead) + missingInRead
		#~ globalTPsum = sum(globalTPlistForARead) + TPsum
		#~ globalFPsum = sum (globalFPlistForARead) + FPsum
		#~ globalRec = globalTPsum / (globalTPsum + globalFNsum) if (globalTPsum + globalFNsum) != 0 else 0
		#~ globalPrec = globalTPsum / (globalTPsum + globalFPsum) if (globalTPsum + globalFPsum) != 0 else 0
		#~ globalPrecision.append(globalPrec)
		#~ globalRecall.append(globalRec)
	GCRateRef.append(GCRateRefRead)
	GCRateCorr.append(GCRateCorrRead)
	return recall, precision,  corBasesRate, uncorCorBasesRate, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, totalCorBases, totalUncorBases


def gapsAndExtensions(reference, corrected, uncorrected, gapsPositions, isExtended, isTrimmed, extendedBasesCount, missingSize):
	refGapsLeft = nbLeftGaps(reference)
	uncoGapsLeft = nbLeftGaps(uncorrected)
	gapsLeft = min(refGapsLeft, uncoGapsLeft)
	if (gapsLeft >= THRESH):
		gapsPositions.extend(i for i in range(gapsLeft))
		if gapsLeft >= THRESH2:
			isExtended = True
			extendedBasesCount.append(gapsLeft - corrected[:gapsLeft].count('.'))
	refGapsRight = nbRightGaps(reference)
	uncoGapsRight = nbRightGaps(uncorrected)
	gapsRight = min(refGapsRight, uncoGapsRight)
	if (gapsRight >= THRESH):
		gapsPositions.extend(i for i in range(len(reference) - 1, len(reference) - gapsRight, - 1))
		if gapsRight >= THRESH2:
			isExtended = True
			extendedBasesCount.append(gapsRight - corrected[len(reference) - gapsRight + 1:].count('.'))
	stretches = findGapStretches(corrected, reference, gapsPositions)
	totalGaps =  gapsLeft + gapsRight
	for s in stretches:
		#~ print(s, stretches[s])
		missingSize += stretches[s] - s - (reference[s: stretches[s]+1].count('.'))
	missingSize -= totalGaps
	if missingSize < 0:
		missingSize = 0
	if missingSize > THRESH:
		isTrimmed = True
	#~ print('missingSize', missingSize)
	return gapsPositions, isExtended , extendedBasesCount, missingSize, stretches, isTrimmed, totalGaps

#~ def nucleotideMetrics(reference, corrected, uncorrected, correctedPositionsRead, existingCorrectedPositionsInThisRead, reportedThreshold, ratioHomopolymers, gapsPositions, indelsubsCorr,  corBasesForARead, uncorBasesForARead, FPlistForARead, TPlistForARead, FNlistForARead, globalFPlistForARead, globalTPlistForARead, globalFNlistForARead, allLenCorrected):
def nucleotideMetrics(reference, corrected, uncorrected,  existingCorrectedPositionsInThisRead, reportedThreshold, ratioHomopolymers, gapsPositions, indelsubsCorr,  corBasesForARead, uncorBasesForARead, uncorCorBasesForARead, uncorUncorBasesForARead, FPlistForARead, TPlistForARead, FNlistForARead,  allLenCorrected, allLenUncorrected):
	FP, TP, FN, corBases, uncorBases, uncorCorBases, uncorUncorBases, GCRateRefRead, GCRateCorrRead, insU, deleU, subsU, insC, deleC, subsC, ratioHomopolymers = getTPFNFP(reference, corrected, uncorrected, existingCorrectedPositionsInThisRead, reportedThreshold, ratioHomopolymers, gapsPositions)
	indelsubsCorr[0] += insC
	indelsubsCorr[1] += deleC
	indelsubsCorr[2] += subsC
	#~ indelsubsUncorr[0] += insU
	#~ indelsubsUncorr[1] += deleU
	#~ indelsubsUncorr[2] += subsU
	corBasesForARead.append(corBases)
	uncorBasesForARead.append(uncorBases)
	uncorCorBasesForARead.append(uncorCorBases)
	uncorUncorBasesForARead.append(uncorUncorBases)
	FPlistForARead.append(FP)
	TPlistForARead.append(TP)
	FNlistForARead.append(FN)
	allLenCorrected.append(getLen(corrected))
	allLenUncorrected.append(getLen(uncorrected))
	return indelsubsCorr,  corBasesForARead, uncorBasesForARead, uncorCorBasesForARead, uncorUncorBasesForARead, FPlistForARead, TPlistForARead, FNlistForARead, allLenCorrected, allLenUncorrected, GCRateRefRead, GCRateCorrRead, insU, deleU, subsU

def computeMetrics(fileName, outPerReadMetrics, correctedFileName, reportedThreshold, clipsNb, readsToSplit):
	msa = open(fileName, 'r')
	nbLines = 0
	lines = msa.readlines()
	nbReadsToDivide = 0
	readNo = 0
	#~ correctedReadsList = getCorrectedReads(correctedFileName)
	#~ upperCasePositions = getUpperCasePositions(correctedReadsList, lines)
	countReadSplit = 0
	countReadExtended = 0
	countReadTrimmed = 0
	extendedBasesCount = []
	missingSize = []
	indelsubsCorr = [0,0,0] #ins del subs
	indelsubsUncorr = [0,0,0] #ins del subs
	allLenCorrected = []
	allLenUncorrected = []
	precision = []
	recall = []
	corBasesRate = []
	uncorCorBasesRate = []
	totalCorBases = 0
	totalUncorBases = 0
	uncorTotalCorBases = 0
	uncorTOtalUncorBases = 0
	totalUncorCorBases = 0
	totalUncorUncorBases = 0
	GCRateRef = []
	GCRateCorr = []
	while nbLines < len(lines):
		if not ">" in lines[nbLines]:
			nbFragments = readsToSplit[headerNo]
			isExtended = False
			isSplit = False
			isTrimmed = False
			gapsPositions = []
			FPlistForARead = []
			TPlistForARead = []
			FNlistForARead = []
			corBasesForARead = []
			uncorBasesForARead = []
			uncorCorBasesForARead = []
			uncorUncorBasesForARead = []
			ratioHomopolymers = []
			missingInRead = 0
			GCRateRefRead = 0
			GCRateCorrRead = 0
			if nbFragments > 1: #split read
				countReadSplit += 1
				splits = 1
				#~ realNotMissing = 0
				realNotMissing = []
				tmpindelsubsUncorr = [[],[],[]]
				while splits <= nbFragments:
					#~ print(splits, nbLines)
					reference = lines[nbLines].rstrip() # get msa for ref
					nbLines += 2
					corrected =  lines[nbLines].rstrip() # msa for uncorrected
					nbLines += 2
					uncorrected = lines[nbLines].rstrip() # msa for corrected
					nbLines += 1
					if len(reference) > 10:
						# upperCasePositions = getUpperCasePositions(correctedFileName, headerNo, corrected)
						# gaps and extensions
						gapsPositions, isExtended , extendedBasesCount, missingInRead, stretches, isTrimmed, totalGaps = gapsAndExtensions(reference, corrected, uncorrected, gapsPositions, isExtended, isTrimmed, extendedBasesCount, missingInRead)
						## zones where the corrected read does not exist / where the correction is not done
						# correctedPositionsRead, existingCorrectedPositionsInThisRead, clips = getCorrectedPositions(stretches, corrected, readNo, upperCasePositions, reference,  clipsNb, header, gapsPositions)
						existingCorrectedPositionsInThisRead, clips = getCorrectedPositions(stretches, corrected, readNo, reference,  clipsNb, header, gapsPositions)
						## indels, subs, TP, FP, FN...
						# indelsubsCorr,  corBasesForARead, uncorBasesForARead, FPlistForARead, TPlistForARead, FNlistForARead, globalFPlistForARead, globalTPlistForARead, globalFNlistForARead, allLenCorrected, GCRateRefRead, GCRateCorrRead,  insU, deleU, subsU  = nucleotideMetrics(reference, corrected, uncorrected, correctedPositionsRead, existingCorrectedPositionsInThisRead, reportedThreshold, ratioHomopolymers, gapsPositions, indelsubsCorr,  corBasesForARead, uncorBasesForARead, FPlistForARead, TPlistForARead, FNlistForARead, globalFPlistForARead, globalTPlistForARead, globalFNlistForARead, allLenCorrected)
						indelsubsCorr,  corBasesForARead, uncorBasesForARead, uncorCorBasesForARead, uncorUncorBasesForARead, FPlistForARead, TPlistForARead, FNlistForARead,  allLenCorrected, allLenUncorrected, GCRateRefRead, GCRateCorrRead,  insU, deleU, subsU  = nucleotideMetrics(reference, corrected, uncorrected,  existingCorrectedPositionsInThisRead, reportedThreshold, ratioHomopolymers, gapsPositions, indelsubsCorr,  corBasesForARead, uncorBasesForARead, uncorCorBasesForARead, uncorUncorBasesForARead, FPlistForARead, TPlistForARead, FNlistForARead, allLenCorrected, allLenUncorrected)
						for pos,cor in enumerate(existingCorrectedPositionsInThisRead):
							if cor:
								realNotMissing.append(pos)
						#~ print(len(corrected), existingCorrectedPositionsInThisRead.count(True), insU, deleU, subsU)
						indelsubsUncorr[0] += insU
						indelsubsUncorr[1] += deleU
						indelsubsUncorr[2] += subsU
						if splits == readsToSplit[headerNo]:
							missingInRead = 0
							for pos in range(len(reference)):
								if pos not in realNotMissing and reference[pos] != '.':
									missingInRead += 1
							maxim = 0
							#~ globalFNlistForARead.append(missingInRead) #add final missed length because of split
							recall, precision,  corBasesRate, uncorCorBasesRate, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, totalCorBases, totalUncorBases = outputMetrics(recall, precision,corBasesRate, uncorCorBasesRate, missingInRead, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, FPlistForARead, TPlistForARead, FNlistForARead, corBasesForARead, uncorBasesForARead, uncorCorBasesForARead, uncorUncorBasesForARead, totalCorBases, totalUncorBases, uncorTotalCorBases, uncorTOtalUncorBases,  GCRateRefRead, GCRateCorrRead)
							if isExtended:
								countReadExtended += 1
							nbReadsToDivide += 1
						if nbLines < len(lines) :
							headerNo = lines[nbLines].split(">")[1].split(" ")[0]
							header = lines[nbLines].split(">")[1].rstrip()
							nbLines += 1
					else:
						if nbLines < len(lines) :
							headerNo = lines[nbLines].split(">")[1].split(" ")[0]
							header = lines[nbLines].split(">")[1].rstrip()
							nbLines += 1
					readNo += 1
					splits += 1

			else: # not split
				
				reference = lines[nbLines].rstrip() # get msa for ref
				nbLines += 2
				corrected =  lines[nbLines].rstrip() # msa for corrected
				nbLines += 2
				uncorrected = lines[nbLines].rstrip()
				nbLines += 1
				if len(reference) > 10:
					# upperCasePositions = getUpperCasePositions(correctedFileName, headerNo, corrected)
					# gaps and extensions
					gapsPositions, isExtended , extendedBasesCount, missingInRead, stretches, isTrimmed, totalGaps = gapsAndExtensions(reference, corrected, uncorrected, gapsPositions, isExtended, isTrimmed, extendedBasesCount, missingInRead)
					## zones where the corrected read does not exist / where the correction is not done
					# correctedPositionsRead, existingCorrectedPositionsInThisRead, clips = getCorrectedPositions(stretches, corrected, readNo, upperCasePositions, reference,  clipsNb, header, gapsPositions)
					existingCorrectedPositionsInThisRead, clips = getCorrectedPositions(stretches, corrected, readNo, reference,  clipsNb, header, gapsPositions)
					## indels, subs, TP, FP, FN...
					# indelsubsCorr,  corBasesForARead, uncorBasesForARead, FPlistForARead, TPlistForARead, FNlistForARead, globalFPlistForARead, globalTPlistForARead, globalFNlistForARead, allLenCorrected, GCRateRefRead, GCRateCorrRead, insU, deleU, subsU = nucleotideMetrics(reference, corrected, uncorrected, correctedPositionsRead, existingCorrectedPositionsInThisRead, reportedThreshold, ratioHomopolymers, gapsPositions, indelsubsCorr, corBasesForARead, uncorBasesForARead, FPlistForARead, TPlistForARead, FNlistForARead, globalFPlistForARead, globalTPlistForARead, globalFNlistForARead, allLenCorrected)
					indelsubsCorr,  corBasesForARead, uncorBasesForARead, uncorCorBasesForARead, uncorUncorBasesForARead, FPlistForARead, TPlistForARead, FNlistForARead,  allLenCorrected, allLenUncorrected, GCRateRefRead, GCRateCorrRead, insU, deleU, subsU = nucleotideMetrics(reference, corrected, uncorrected,  existingCorrectedPositionsInThisRead, reportedThreshold, ratioHomopolymers, gapsPositions, indelsubsCorr, corBasesForARead, uncorBasesForARead, uncorCorBasesForARead, uncorUncorBasesForARead, FPlistForARead, TPlistForARead, FNlistForARead,  allLenCorrected, allLenUncorrected)
					indelsubsUncorr[0] += insU
					indelsubsUncorr[1] += deleU
					indelsubsUncorr[2] += subsU
					recall, precision, corBasesRate, uncorCorBasesRate, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, totalCorBases, totalUncorBases= outputMetrics(recall, precision,corBasesRate, uncorCorBasesRate, missingInRead, missingSize, GCRateRef, GCRateCorr, outPerReadMetrics, FPlistForARead, TPlistForARead, FNlistForARead, corBasesForARead, uncorBasesForARead, uncorCorBasesForARead, uncorUncorBasesForARead, totalCorBases, totalUncorBases, uncorTotalCorBases, uncorTOtalUncorBases, GCRateRefRead, GCRateCorrRead)
					if isExtended:
						countReadExtended += 1
					if isTrimmed:
						countReadTrimmed += 1
					nbReadsToDivide += 1
					readNo += 1
					if nbLines < len(lines):
						headerNo = lines[nbLines].split(">")[1].split(" ")[0]
						header = lines[nbLines].split(">")[1].rstrip()
						nbLines += 1
				else:
					readNo += 1
					if nbLines < len(lines):
						headerNo = lines[nbLines].split(">")[1].split(" ")[0]
						header = lines[nbLines].split(">")[1].rstrip()
						nbLines += 1
		else:
			headerNo = lines[nbLines].split(">")[1].split(" ")[0]
			header = lines[nbLines].split(">")[1].rstrip()
			nbLines += 1

	# compute global metrics
	GCRateRef = round(sum(GCRateRef) / len(GCRateRef),3)
	GCRateCorr = round(sum(GCRateCorr) / len(GCRateCorr),3)
	recall = sum(recall)*1.0 / nbReadsToDivide if nbReadsToDivide != 0 else 0
	precision = sum(precision)*1.0 / nbReadsToDivide if nbReadsToDivide != 0 else 0
	#~ globalRecall = sum(globalRecall)*1.0 / nbReadsToDivide if nbReadsToDivide != 0 else 0
	#~ globalPrecision = sum(globalPrecision)*1.0 / nbReadsToDivide if nbReadsToDivide != 0 else 0
	corBasesRate = sum(corBasesRate)*1.0 / nbReadsToDivide if nbReadsToDivide != 0 else 0
	uncorCorBasesRate = sum(uncorCorBasesRate)*1.0 / nbReadsToDivide if nbReadsToDivide != 0 else 0
	throughput = sum(allLenCorrected)
	uncorThroughput = sum (allLenUncorrected)
	errorRate = 1 - (totalCorBases / (totalCorBases + totalUncorBases))
	uncorErrorRate = 1 - (totalUncorBases / (totalCorBases + totalUncorBases))
	if len(ratioHomopolymers) > 1:
		meanRatioHomopolymers = statistics.mean(ratioHomopolymers)
	else:
		meanRatioHomopolymers = 1
	return nbReadsToDivide, throughput, uncorThroughput, precision, recall, corBasesRate, errorRate, uncorCorBasesRate, uncorErrorRate, missingSize,  GCRateRef, GCRateCorr, indelsubsUncorr, indelsubsCorr, meanRatioHomopolymers, allLenCorrected, countReadSplit, countReadTrimmed, countReadExtended, extendedBasesCount


# get the position of nt in uppercase to compute recall and precision only at these positions
#~ def getUpperCasePositions(correctedReadsList, lines):
	#~ upperCasePositions = [] # positions to take into account in the msa
	#~ nbLines = 3 # starting at the first corrected sequence
	#~ headerNo = lines[0].split(">")[1].split(" ")[0]
	#~ headerPrev = ""
	
	#~ while nbLines < len(lines):
		#~ headerNo = lines[nbLines - 1].split(">")[1].split(" ")[0]
		#~ correctedMsa = lines[nbLines].rstrip()
		#~ #TOVERIFY
		#~ if headerNo == headerPrev:
			#~ index += 1
		#~ else:
			#~ index = 0
		#~ if headerNo in correctedReadsList.keys():
			#~ correctedReadSequence = correctedReadsList[headerNo][index]
			#~ headerPrev = headerNo
			#~ posiNt = 0
			#~ posiNtSeq = 0
			#~ inUpper = False
			#~ upperCasePositions.append([])
			#~ while posiNt < len(correctedMsa):
				#~ if posiNtSeq >= len(correctedReadSequence):
					#~ upperCasePositions[-1].append(False)
					#~ posiNt += 1
				#~ else:
					#~ nt = correctedMsa[posiNt]
					#~ ntSeq = correctedReadSequence[posiNtSeq]
					#~ if not ntSeq.islower():
						#~ if nt != ".":
							#~ inUpper = True
						#~ else:
							#~ inUpper = False
					#~ else:
						#~ inUpper = False
					
							
					#~ if inUpper:
						#~ upperCasePositions[-1].append(True)
					#~ else:
						#~ upperCasePositions[-1].append(False)
					#~ if nt != ".":
						#~ posiNtSeq += 1
					#~ posiNt += 1
			#~ nbLines += 6
		#~ else:
			#~ upperCasePositions[-1] = [False] * len(correctedMsa)
	#~ return upperCasePositions


def getUpperCasePositions(correctedReadsFile, header, correctedMsa):
	correctedReadSequence = utils.getCorrectedSequence(correctedReadsFile, header)
	upperCasePositions = [] # positions to take into account in the msa
	posiNt = 0
	posiNtSeq = 0
	inUpper = False
	while posiNt < len(correctedMsa):
		if posiNtSeq >= len(correctedReadSequence):
			upperCasePositions.append(False)
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
				upperCasePositions.append(True)
			else:
				upperCasePositions.append(False)
			if nt != ".":
				posiNtSeq += 1
			posiNt += 1
	return upperCasePositions

# add to the uppercase positions the positions where there is no stretch of "." , i.e. all positions where recall and precision are actually computed

# def getCorrectedPositions(stretches, corrected, readNo, upperCasePositions, reference, clipsNb, header,  gapsPositions):
def getCorrectedPositions(stretches, corrected, readNo, reference, clipsNb, header,  gapsPositions):
	msaLineLen = len(corrected)
	#~ correctedPositions = copy.copy(upperCasePositions) #all positions in upper case in the corrected read
	existingCorrectedPositions = [True] * msaLineLen
	leftClipping = 0
	rightClipping = None
	if header in clipsNb.keys():
		leftClipping = clipsNb[header][0]
		rightClipping = msaLineLen - clipsNb[header][1]

	i = 0
	j = 0
	lenClip = 0
	while j < leftClipping:
		if corrected[i] != ".":
			j += 1
		existingCorrectedPositions[i] = False
		i += 1
		lenClip += 1
	#~ print(i)

	i = msaLineLen - 1
	j = msaLineLen - 1
	if rightClipping is not None:
		while j >= rightClipping:
			if corrected[i] != ".":
				j -= 1
			existingCorrectedPositions[i] = False
			lenClip += 1
			i -= 1

	#~ for pos in forNotExisting:
		#~ existingCorrectedPositions[pos] = False

	positionsToRemove = list()
	if len(stretches.keys()) > 0:  # split read (or trimmed)
		for pos in stretches.keys(): 
			positionsToRemove.append([pos, stretches[pos]]) #interval(s)) in which the corrected read sequence does not exist
		for interv in positionsToRemove:
			for i in range(interv[0], interv[1] +1):
				#~ print(interv[0], interv[1])
				existingCorrectedPositions[i] = False # remove regions where there is no corrected sequence  # remove regions where there is no corrected sequence (split/trimmed) from corrected regions
				#~ correctedPositions[i] = False
	for i in gapsPositions:
		existingCorrectedPositions[i] = False
		#~ correctedPositions[i] = False
	#~ return correctedPositions, existingCorrectedPositions, lenClip
	return existingCorrectedPositions, lenClip
