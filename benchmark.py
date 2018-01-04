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
#~ import readAndSortFiles
#~ import computeStats
import alignment
import computeStats
import readAndSortFiles



# count the number of reads in a file
def getFileReadNumber(fileName):
	cmdGrep = """grep ">" -c """ + fileName
	val = subprocess.check_output(['bash','-c', cmdGrep])
	return int(val.decode('ascii'))



def main():
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	installDirectory = os.path.dirname(os.path.realpath(__file__))
	# Manage command line arguments
	parser = argparse.ArgumentParser(description="Benchmark for quality assessment of long reads correctors.")
	# Define allowed options
	corrected = ""
	uncorrected = ""
	reference = ""
	parser = argparse.ArgumentParser()
	parser.add_argument('-threads', nargs='?', type=int, action="store", dest="threads", help="Number of threads", default=2)
	parser.add_argument('-corrected', nargs='?', type=str, action="store", dest="corrected", help="Fasta file with corrected reads (each read sequence on one line)")
	parser.add_argument('-uncorrected', nargs='?', type=str,  action="store", dest="uncorrected",  help="Fasta file with uncorrected reads (each read sequence on one line)")
	parser.add_argument('-reference', nargs='?', type=str,  action="store", dest="reference",  help="Fasta file with reference read sequences (each read sequence on one line)")
	parser.add_argument('-tool', nargs='?', type=str,  action="store", dest="soft",  help="Corrector used (lowercase, in this list: lorma, mecat, pbdagcon, daccord). If no corrector name is provided, make sure the read's headers are correctly formatted (i.e. they correspond to those of uncorrected and reference files)")
	# get options for this run
	args = parser.parse_args()
	if (len(sys.argv) <= 1):
		parser.print_help()
		return 0
	
	
	corrected = args.corrected
	uncorrected = args.uncorrected
	reference = args.reference #todo : won't be necessary when a real simulator will be used
	soft = None
	if args.soft is not None:
		if args.soft == "lorma" or args.soft == "mecat" or args.soft == "pbdagcon" or args.soft == "daccord":
			soft = args.soft
	size =  getFileReadNumber(corrected)
	readAndSortFiles.processReadsForAlignment(soft, reference, uncorrected, corrected, size, soft)
	alignment.getPOA(corrected, reference, uncorrected, args.threads, installDirectory, soft)
	computeStats.outputRecallPrecision(corrected, 0, 0, soft)


if __name__ == '__main__':
	main()

