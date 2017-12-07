#!/usr/bin/env python3

"""*****************************************************************************
 *   Authors: Camille Marchet  Pierre Morisse Lolita Lecompte Antoine Limasset
 *   Contact: camille.marchet@irisa.fr, IRISA/Univ Rennes/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
 *   Source: https://github.com/Malfoy/BGREAT
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
	#~ return sortedList
	outfile = open(outfileName, 'w')
	for record in sortedList:
		outfile.write(">" + record.description+"\n")
		outfile.write(str(record.seq)+"\n")
	#~ print(sortedList)

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
	outfile = "corrected_by_lordec.fa"
	logFile = open("lordec.log", 'w')
	cmdLordec = "lordec-correct -i simulatedReads.fa -2 simulatedReads_short.fa -o " + outfile + " -k 21 -s 2 -T " + str(threads)
	p = subprocessLauncher(cmdLordec, logFile, logFile)
	cmdRm = "rm *.h5"
	subprocess.check_output(['bash','-c', cmdRm])
	return outfile


def colormap(threads):
	outfile = "corrected_by_colormap.fa"
	logFile = open("colormap.log", 'w')
	cmdColorMap = "runCorr.sh simulatedReads.fa simulatedReads_short.fa colorMapDir pre " + str(threads)
	cmdColorMap = "runOEA.sh colorMapDir/pre_sp.fasta simulatedReads_short.fa colorMapDir pre " + str(threads)
	p = subprocessLauncher(cmdColorMap, logFile, logFile)
	cmdCp = "cp colorMapDir/pre_oea.fasta " + outfile
	subprocess.check_output(['bash','-c', cmdCp])
	return outfile


def lorma(threads):
	cmdcp = "cp simulatedReads.fa" + suffix + ".fa copy.fa"
	subprocess.check_output(['bash','-c', cmdcp])
	cmdLorma = "lorma.sh copy.fa -s -threads " + str(threads)
	p = subprocessLauncher(cmdLorma)
	try:
		cmdmv = "mv trim.fasta corrected_by_LoRMA_trim.fa"
		subprocess.check_output(['bash','-c', cmdmv])
	except subprocess.CalledProcessError:
		pass
	try:
		cmdmv = "mv discarded.fasta discarded_by_LoRMA_discarded.fa"
		subprocess.check_output(['bash','-c', cmdmv])
	except subprocess.CalledProcessError:
		pass
	try:
		cmdmv = "mv final.fasta corrected_by_LoRMA.fa"
		subprocess.check_output(['bash','-c', cmdmv])
	except subprocess.CalledProcessError:
		cmdnull = "> corrected_by_LoRMA" + suffix + ".fa"
		subprocess.check_output(['bash','-c', cmdnull])

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



def computeMetrics(fileName, outfile):
	msa = open(fileName, 'r')
	sumFP = []
	sumFN = []
	sumTP = []
	nbLines = 0
	lines = msa.readlines()
	readNo = 0
	while nbLines < len(lines) - 3:
		toW = ""
		if not ">" in lines[nbLines]:
			reference = lines[nbLines].rstrip()
			nbLines += 2
			uncorrected =  lines[nbLines].rstrip()
			nbLines += 2
			corrected = lines[nbLines].rstrip()
			nbLines += 2
			FN = 0 #code M
			FP = 0 #code !
			TP = 0 #code *
			for ntRef, ntUnco, ntResult in zip(reference, uncorrected, corrected):
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
	return precision, recall







def main():
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	# Manage command line arguments
	parser = argparse.ArgumentParser(description="Benchmark for quality assessment of long reads correctors.")
	# Define allowed options
	threads = 2
	parser = argparse.ArgumentParser()
	parser.add_argument('-genome', nargs='?', type=str, action="store", dest="genomeRef", help="Reference genome file for simulation (sequence on one line)")
	parser.add_argument('-read_length', nargs='?', type=int, action="store", dest="readLen", help="Simulated read length", default=10000)
	parser.add_argument('-coverage', nargs='?', type=int, action="store", dest="coverage", help="Simulation coverage (example: 10 for 10X)", default=10)
	parser.add_argument('-error_rate', nargs='?', type=float, action="store", dest="errorRate", help="Error rate (example: 0.1 for 10 percent)", default = 0.1)
	parser.add_argument('-threads', nargs='?', type=int, action="store", dest="threads", help="Number of threads", default = 2)
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

	if args.corrected is None and args.uncorrected is None and args.reference is None:
		# simulate data
		cmdSimul = "./bin/simulator " + args.genomeRef +  " " + str(args.readLen) + " " + str(args.coverage) + " " + str(args.errorRate) + " simulatedReads "
		uncorrected = "simulatedReads.fa"
		reference = "p.simulatedReads.fa"
		subprocessLauncher(cmdSimul)
		# launch correctors
		# we assume binaries are in PATH
		for soft in ["lordec", "colormap", "mecat", "lorma"]:
			if soft == "lordec":
				beg = time.time()
				corrected = lordec(threads)
				end = time.time()
			#~ if soft == "mecat":
				#~ beg = time.time()
				#~ corrected = mecat()
				#~ end = time.time()
			if soft == "colormap":
				beg = time.time()
				corrected_tmp = colormap(threads)
				end = time.time()
				corrected = "corrected_sorted_by_colormap.fa"
				readAndSortFasta(corrected_tmp, corrected)
			outProfile = open(soft + "_msa_profile.txt", 'w')
			cmdPOA = "./bin/poa -corrected_reads_fasta " + corrected + " -reference_reads_fasta " + reference + " -uncorrected_reads_fasta " + uncorrected + "-threads " + str(threads)
			subprocessLauncher(cmdPOA)
			# gets precision and recall from MSA of 3 versions of reads
			cmdMv = "mv default_output_msa.fasta msa_" + soft + ".fa"
			subprocess.check_output(['bash','-c', cmdMv])
			precision, recall = computeMetrics("msa_" + soft + ".fa", outProfile)
			outProfile.write("\n***********SUMMARY***********\n")
			outProfile.write(soft + ": Recall " + str(round(recall,2)) + " Precision " + str(round(precision,2)) + "\n")
			print(soft + ": Recall ", round(recall,2), "Precision ", round(precision,2))
			outProfile.write("Run in {0} seconds.".format(str(round(end-beg, 2)))+"\n") #runtime of the tool
			print(soft + " ran in {0} seconds.".format(str(round(end-beg, 2)))) #runtime of the tool
		
	else: # else directly use data provided and skip simulation
		corrected = args.corrected
		uncorrected = args.uncorrected
		reference = args.reference
		outProfile = open("msa_profile.txt", 'w')
		cmdPOA = "./bin/poa -corrected_reads_fasta " + corrected + " -reference_reads_fasta " + reference + " -uncorrected_reads_fasta " + uncorrected + " -threads " + str(threads)
		subprocessLauncher(cmdPOA)
		# gets precision and recall from MSA of 3 versions of reads
		cmdMv = "mv default_output_msa.fasta msa.fa"
		subprocess.check_output(['bash','-c', cmdMv])
		precision, recall = computeMetrics("msa.fa", outProfile)
		outProfile.write("\n***********SUMMARY***********\n")
		outProfile.write(": Recall " + str(round(recall,2)) + " Precision " + str(round(precision,2)) + "\n")
		print( "Recall ", round(recall,2), "Precision ", round(precision,2))
		#~ outProfile.write("Run in {0} seconds.".format(str(round(end-beg, 2)))+"\n") #runtime of the tool
		#~ print(" ran in {0} seconds.".format(str(round(end-beg, 2)))) #runtime of the tool
		# launch poa graph for MSA: prerequisite = all the sequences file have the same size and sequences come in the same order
		


if __name__ == '__main__':
	main()

