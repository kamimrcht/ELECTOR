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

#Assemble the reads contained in the file readsF, with nbThreads threads
def runAssembly(readsF, nbThreads):
	#align and assemble the reads
	reFile = basename(os.path.splitext(readsF)[0])
	cmdAl = "minimap -Sw5 -L100 -m0 -t " + nbThreads + " "  + readsF + " " + readsF
	outErr = open("/dev/null", 'w')
	alFile = reFile + ".paf"
	outAl = open(alFile, 'w')
	subprocessLauncher(cmdAl, outAl, outErr)
	outAl.close()
	cmdAs = "miniasm -f " + readsF + " " + alFile
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
	cmdRm = "rm " + alFile + " " + asFile
	subprocessLauncher(cmdRm)

#Compute the number of contigs and the NG50 of the assembly contained
#in the file contigs, based on the provided genome size genSize.
def computeContigsNbAndNG50(contigs, genSize):
	f = open(contigs)
	sizes = []
	line = f.readline()
	while line != '':
		line = f.readline()[:-1]
		sizes.append(len(line))
		line = f.readline()
	f.close()
	sortedSizes = [f for f in sorted(sizes, key=lambda x : int(x), reverse=True)]
	ng50 = 0
	i = 0
	while i < len(sortedSizes) and ng50 < 1 / 2 * genSize:
		ng50 = ng50 + sortedSizes[i]
		i = i + 1
	return [len(sortedSizes),sortedSizes[i-1]]

def computeCovAndId(contigs, reference):
	cmdAl = "dnadiff " + reference + " " + contigs
	out = open("/dev/null", 'w')
	subprocessLauncher(cmdAl, out, out)
	out.close()
	f = open("out.report")
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	cov = line.split("(")[1].split(")")[0]
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	id = line.split(" ")[-1:][0][:-1]
	cmdRm = "rm out.1coords out.1delta out.delta out.mcoords out.mdelta out.qdiff out.rdiff out.report out.snps"
	subprocessLauncher(cmdRm)
	f.close()
	return [cov,id]

#runAssembly(sys.argv[1], sys.argv[3])
nbContigsNG50 = computeContigsNbAndNG50(basename(os.path.splitext(sys.argv[1])[0]) + ".contigs.fa", getTotalLength(sys.argv[2]))
nbContigs = nbContigsNG50[0]
NG50 = nbContigsNG50[1]
covId = computeCovAndId(sys.argv[1], sys.argv[2])
cov = covId[0]
id = covId[1]
print("Number of contigs : " + str(nbContigs))
print("NG50 : " + str(NG50))
print("Genome coverage : " + cov)
print("Identity : " + str(id) + "%")
