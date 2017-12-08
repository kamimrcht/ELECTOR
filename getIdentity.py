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
		line = f.readline()
		if line[0] != '>':
			totalLength = totalLength + len(line)
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
			nbs = [int(i) for i in (re.findall('\d+', q))]
			out.write(str(sum(nbs) / l * 100) + '\n')
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


readsBaseName = sys.argv[2].split(".")[0]
cmdAlign = "bwa mem -t " + sys.argv[3] + " " + sys.argv[1] + " " + sys.argv[2]
outSam = open(readsBaseName + ".sam", 'w')
outDevNull = open("/dev/null", 'w')
subprocessLauncher(cmdAlign, outSam, outDevNull)
outSam.close()
computeIdentity(readsBaseName + ".sam", readsBaseName + ".id")
avId = averageIdentity(readsBaseName + ".id")
cmdConvertToBam = "samtools view -Sb " + readsBaseName + ".sam"
outBam = open(readsBaseName + ".bam", 'w')
cmdSortBam = "samtools sort " + readsBaseName + ".bam"
outSBam = open("sorted_" + readsBaseName + ".bam", 'w')
cmdGetCov = "samtools depth sorted_" + readsBaseName + ".bam"
outCov = open(readsBaseName + ".cov", 'w')
subprocessLauncher(cmdConvertToBam, outBam)
outBam.close()
subprocessLauncher(cmdSortBam, outSBam)
subprocessLauncher(cmdGetCov, outCov)
outCov.close()
refLength = getTotalLength(sys.argv[1])
inCov = open(readsBaseName + ".cov")
coveredBases = sum(1 for line in inCov)
inCov.close()

print("Average identity : ", avId, "%")
print("Genome covered : ", float(coveredBases / refLength * 100), "%")
