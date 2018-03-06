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
import alignment
from multiprocessing import Pool, TimeoutError


installDirectoryGlobal=""


try:
    from subprocess import DEVNULL # py3k
except ImportError:
    import os
    DEVNULL = open(os.devnull, 'wb')



# launch subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
	args = shlex.split(cmd)
	#~ p = subprocess.Popen(args, stdin = argstdin, stdout = argstdout, stderr = argstderr).communicate()
	p = (subprocess.Popen(args, stdin = argstdin, stdout = DEVNULL, stderr = DEVNULL))
	#~ rc=p.communicate()
	#~ print (p.wait())
	return p.wait()

#~ # do the msa
#~ def getPOA(corrected, reference, uncorrected, threads, installDirectory, soft=None):

def f(i):
	cmdPOA = installDirectoryGlobal + "/bin/poa -pir swag"+str(i)+"  -preserve_seqorder -corrected_reads_fasta out3"+str(i)+" -reference_reads_fasta out1"+str(i)+" -uncorrected_reads_fasta out2"+str(i)+" -preserve_seqorder -threads  1 -pathMatrix " + installDirectoryGlobal
	#~ print(cmdPOA)
	subprocessLauncher(cmdPOA)
	return i


def getPOA(corrected, reference, uncorrected, threads, installDirectory, outDir, soft=None):
	oldMode=False
	#~ oldMode=True
	if(oldMode):
		cmdPOA = installDirectory + "/bin/poa -preserve_seqorder -corrected_reads_fasta " + corrected + " -reference_reads_fasta " + reference + " -uncorrected_reads_fasta " + uncorrected + " -threads " + str(threads) + "  -pathMatrix " + installDirectory
		subprocessLauncher(cmdPOA)
		if soft is not None:
			cmdMv = "mv default_output_msa.fasta " + outDir + "/msa_" + soft + ".fa"
		else:
			cmdMv = "mv default_output_msa.fasta " + outDir + "/msa.fa"
		subprocess.check_output(['bash','-c', cmdMv])
	else:
		amount_nuc=100*1000*1000;
		print("- mean that a large amount of nuc has been handled (100,000,000)")
		global installDirectoryGlobal
		installDirectoryGlobal=installDirectory
		position_in_read_file=1

		cmdRM = "rm progress.txt"
		subprocess.call(['bash','-c', cmdRM],stdout=DEVNULL,stderr=DEVNULL)
		while(position_in_read_file!=0):
			cmdSplitter = installDirectory + "/bin/masterSplitter "+ reference +" "+uncorrected+" "+corrected +" out1 out2 out3 15 100 "+str(amount_nuc)
			#~ print(cmdSplitter)
			position_in_read_file=subprocessLauncher(cmdSplitter)
			#~ print(position_in_read_file)
			#~ print("Wait for 100 '-' to be printed")
			with Pool (processes=threads) as pool:
				for i in pool.imap_unordered(f, range(100)):
					continue
					#~ sys.stdout.write('-')
					#~ sys.stdout.flush()
			#~ print()
				#~ print(pool.map(f, range(100)))
			#~ for i in range(0, 100):
				#~ for(j in range(0,threads)):
					#~ cmdPOA = installDirectory + "/bin/poa -pir swag"+str(i)+"  -preserve_seqorder -corrected_reads_fasta out3"+str(i)+" -reference_reads_fasta out1"+str(i)+" -uncorrected_reads_fasta out2"+str(i)+" -preserve_seqorder -threads  1 -pathMatrix " + installDirectory
			for i in range(0, 100):
				cmdMerger = installDirectory + "/bin/Donatello swag"+str(i)+"  merger"
				subprocessLauncher(cmdMerger)
			sys.stdout.write('-')
			sys.stdout.flush()
			cmdRM = "rm out*"
			subprocess.call(['bash','-c', cmdRM])
			cmdRM = "rm swag*"
			subprocess.call(['bash','-c', cmdRM])

		print()
		if soft is not None:
			cmdMv = "mv merger " + outDir + "/msa_" + soft + ".fa"
		else:
			cmdMv = "mv merger " + outDir + "/msa.fa"
		subprocess.call(['bash','-c', cmdMv])
		#~ cmdRM = "rm out*"
		#~ subprocess.call(['bash','-c', cmdRM])
		#~ cmdRM = "rm swag*"
		#~ subprocess.call(['bash','-c', cmdRM])
