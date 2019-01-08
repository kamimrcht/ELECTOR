#!/usr/bin/python3

import sys
import argparse
import os
import re
import csv
import shlex, subprocess
from subprocess import Popen, PIPE, STDOUT

#Launches subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None, argstdin=None):
	args = shlex.split(cmd)
	p = subprocess.Popen(args, stdin = argstdin, stdout = argstdout, stderr = argstderr).communicate()
	return p

#Returns the total length of the sequences contained in reference.
def getTotalLength(reference):
	totalLength = 0
	f = open(reference)
	line = f.readline()
	while line != '':
		if line[0] != '>':
			totalLength = totalLength + len(line[:-1])
		line = f.readline()
	f.close()
	return totalLength


#Computes the idetity of each alignement in the file alignements.
#Stores the results in the file ids.
def computeIdentity(alignments, ids):
	f = open(alignments)
	out = open(ids, 'w')
	line = f.readline()
	#Skip headers
	while line[0] == "@":
		line = f.readline()
	while line != '':
		t = line.split("\t")
		#Compute identity only for full alignments
		if t[1] == "0" or t[1] == "16":
			pos = t[3]
			l = len(t[9])
			q = t[12]
			cigar = t[5]
			nbs = [int(i) for i in (re.findall('\d+', q))]
			dels = sum([int(i.split("D")[0]) for i in (re.findall('\d+D', cigar))])
			clips = sum([int(i.split("S")[0]) for i in (re.findall('\d+S', cigar))])
			out.write(str(sum(nbs) / (l+dels-clips) * 100) + '\n')
		line = f.readline()
	f.close()
	out.close()

#Computes the average identity of the alignements in the file alignments
def averageIdentity(alignments):
	f = open(alignments)
	nbReads = 0
	avId = 0
	s = f.readline()
	while s != '':
		nbReads = nbReads + 1
		avId = avId + float(s)
		s = f.readline()
	f.close()
	return avId / nbReads

#Compute the genome coverage of the alignments
def computeCoverage(readsBaseName, reference):
	cmdConvertToBam = "./samtools/samtools view -Sb " + readsBaseName + ".sam"
	outBam = open(readsBaseName + ".bam", 'w')
	cmdSortBam = "./samtools/samtools sort " + readsBaseName + ".bam"
	outSBam = open(readsBaseName + "_sorted.bam", 'w')
	cmdGetCov = "./samtools/samtools depth " + readsBaseName + "_sorted.bam"
	outCov = open(readsBaseName + ".cov", 'w')
	subprocessLauncher(cmdConvertToBam, outBam)
	outBam.close()
	subprocessLauncher(cmdSortBam, outSBam)
	outSBam.close()
	subprocessLauncher(cmdGetCov, outCov)
	outCov.close()
	refLength = getTotalLength(reference)
	inCov = open(readsBaseName + ".cov")
	coveredBases = sum(1 for line in inCov)
	inCov.close()
	cov = float(coveredBases / refLength * 100)

def generateResults(reads, reference, threads, logFile):
	threads = str(threads)

	readsBaseName = os.path.splitext(reads)[0]
	cmdMakeIndex = "./bwa/bwa index " + reference
	cmdAlign = "./bwa/bwa mem -t " + threads + " " + reference + " " + reads
	outSam = open(readsBaseName + ".sam", 'w')
	outErr = open("/dev/null", 'w')
	subprocessLauncher(cmdMakeIndex, outErr, outErr)
	subprocessLauncher(cmdAlign, outSam, outErr)
	outSam.close()
	outErr.close()
	computeIdentity(readsBaseName + ".sam", readsBaseName + ".id")
	avId = averageIdentity(readsBaseName + ".id")
	cov = computeCoverage(readsBaseName, reference)

	print("Average identity : " + str(round(avId, 3)) + "%")
	print("Genome covered : " + str(round(cov, 3)) + "%")
	logFile.write("Average identity : " + str(round(avId, 3)) + "%\nGenome covered : " + str(round(cov, 3)) + "%\n")
	return str(avId), str(cov)
