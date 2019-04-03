#!/usr/bin/env python3

import sys
import argparse
import os
import re
import csv
import shlex, subprocess
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

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

#Assemble the reads contained in the file readsF, with nbThreads threads,
#and return the number of contigs
def runAssembly(readsF, nbThreads):
	#align and assemble the reads
	reFile = (os.path.splitext(readsF)[0])
	cmdAl = "./minimap2/minimap2 -x ava-ont -t " + nbThreads + " "  + readsF + " " + readsF
	outErr = open("/dev/null", 'w')
	alFile = reFile + ".paf"
	outAl = open(alFile, 'w')
	subprocessLauncher(cmdAl, outAl, outErr)
	outAl.close()
	cmdAs = "./miniasm/miniasm -f " + readsF + " " + alFile
	asFile = reFile + ".gfa"
	outAs = open(asFile, 'w')
	subprocessLauncher(cmdAs, outAs, outErr)
	outErr.close()
	outAs.close()

	#format the output file
	cmdSed = "sed -i '/^S/!d' " + asFile
	subprocessLauncher(cmdSed)
	contigsFile = reFile + ".contigs.fa"
	outContigs = open(contigsFile, 'w')
	inAs = open(asFile)
	nbContigs = 1
	line = inAs.readline()
	while line != '':
		seq = line.split("\t")[2]
		outContigs.write(">contig" + str(nbContigs) + "\n" + seq + "\n")
		nbContigs = nbContigs + 1
		line = inAs.readline()
	inAs.close()
	outContigs.close()
	return nbContigs - 1

#Align the contigs to the reference genome
def alignContigs(contigs, reference, nbThreads):
	cmdAl = "minimap2 -a --MD -t " + nbThreads + " " + reference + " " + contigs + ".contigs.fa"
	outAl = open(contigs + ".contigs.sam", 'w')
	outEr = open("/dev/null", 'w')
	subprocessLauncher(cmdAl, outAl, outEr)
	outAl.close()

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
			#q = t[12]
			q = line.split("MD:Z:")[1]
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

#Compute the number of aligned contigs and the NGA50 and NGA75 of the aligned contigs contained
#in the file contigs, based on the provided genome size genSize.
def computeContigsNbAndNG50(contigs, genSize):
	f = open(contigs)
	sizes = []
	line = f.readline()
	while line != '' and line[0] == "@":
		line = f.readline()
	while line != '':
		fields = line.split("\t")
		if fields[1] == "0" or fields[1] == "16":
			sizes.append(len(fields[9]))
		line = f.readline()
	f.close()
	sortedSizes = [f for f in sorted(sizes, key=lambda x : int(x), reverse=True)]
	nga50s = 0
	nga75s = 0
	i = 0
	while i < len(sortedSizes) and nga50s < 1 / 2 * genSize:
		nga50s = nga50s + sortedSizes[i]
		nga75s = nga75s + sortedSizes[i]
		i = i + 1
	nga50 = i - 1
	while i < len(sortedSizes) and nga75s < 75 / 100 * genSize:
		nga75s = nga75s + sortedSizes[i]
		i = i + 1
	nga75 = i - 1
	if len(sortedSizes) > 0:
		return [len(sortedSizes),sortedSizes[nga50],sortedSizes[nga75]]
	else:
		return [0,0,0]

#Computes the number of breakpoints in the assembly.
def computeNbBreakpoints(file):
	cmd = "./samtools/samtools flagstat " + file + ".contigs.sam"
	out = open(file + ".contigs.fs", 'w')
	subprocessLauncher(cmd, out)
	out.close()
	fs = open(file + ".contigs.fs")
	l = fs.readline()
	l = fs.readline()
	l = fs.readline().split(" ")
	bp = int(l[0]) + int(l[2])
	return bp

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
	return cov

def generateResults(reads, reference, threads, logFile):
	threads = str(threads)

	readsBaseName = os.path.splitext(reads)[0]
	nbContigs = runAssembly(reads, threads)
	alignContigs((os.path.splitext(reads)[0]), reference, threads)
	nbContigsNGs = computeContigsNbAndNG50((os.path.splitext(reads)[0]) + ".contigs.sam", getTotalLength(reference))
	nbAlContigs = nbContigsNGs[0]
	NG50 = nbContigsNGs[1]
	NG75 = nbContigsNGs[2]
	nbBreakpoints = computeNbBreakpoints((os.path.splitext(reads)[0]))
	cov = computeCoverage(readsBaseName + ".contigs", reference)
	computeIdentity(readsBaseName + ".contigs.sam", readsBaseName + ".contigs.id")
	id = averageIdentity(readsBaseName + ".contigs.id")

	print("Number of contigs : " + str(nbContigs))
	print("Number of aligned contigs : " + str(nbAlContigs))
	print("Number of breakpoints : " + str(nbBreakpoints))
	print("NGA50 : " + str(NG50))
	print("NGA75 : " + str(NG75))
	print("Genome covered : " + str(round(cov, 4)) + "%")
	print("Identity : " + str(round(id, 4)) + "%")
	logFile.write("Number of contigs : " + str(nbContigs) + "\n" + "Number of aligned contigs : " + str(nbAlContigs) + "\n" + "Number of breakpoints : " + str(nbBreakpoints) + "\n" + "NGA50 : " + str(NG50) + "\n" + "NGA75 : " + str(NG75) + "\n" + "Genome covered : " + str(round(cov, 4)) + "%\n" + "Identity : " + str(round(id, 4)) + "%\n")
	return str(nbContigs), str(nbAlContigs), str(nbBreakpoints),  str(NG50), str(NG75), str(cov)
