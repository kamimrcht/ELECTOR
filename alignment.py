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


try:
    from subprocess import DEVNULL # py3k
except ImportError:
    import os
    DEVNULL = open(os.devnull, 'wb')



# launch subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
	args = shlex.split(cmd)
	#~ p = subprocess.Popen(args, stdin = argstdin, stdout = argstdout, stderr = argstderr).communicate()
	p = subprocess.Popen(args, stdin = argstdin, stdout = DEVNULL, stderr = argstderr).communicate()
	return p

#~ # do the msa
#~ def getPOA(corrected, reference, uncorrected, threads, installDirectory, soft=None):


def getPOA(corrected, reference, uncorrected, threads, installDirectory, soft=None):
	#~ oldMode=False
	oldMode=True
	if(oldMode):
		cmdPOA = installDirectory + "/bin/poa -preserve_seqorder -corrected_reads_fasta " + corrected + " -reference_reads_fasta " + reference + " -uncorrected_reads_fasta " + uncorrected + " -threads " + str(threads) + "  -pathMatrix " + installDirectory
		subprocessLauncher(cmdPOA)
		if soft is not None:
			cmdMv = "mv default_output_msa.fasta msa_" + soft + ".fa"
		else:
			cmdMv = "mv default_output_msa.fasta msa.fa"
		subprocess.check_output(['bash','-c', cmdMv])
	else:
		cmdSplitter = installDirectory + "/bin/masterSplitter "+ reference +" "+uncorrected+" "+corrected +" out1 out2 out3 15 100"
		subprocessLauncher(cmdSplitter)
		for i in range(0, 100):
			print(str(i)+"%")
			cmdPOA = installDirectory + "/bin/poa -corrected_reads_fasta out3"+str(i)+" -reference_reads_fasta out1"+str(i)+" -uncorrected_reads_fasta out2"+str(i)+" -preserve_seqorder -threads  1 -pathMatrix " + installDirectory
			cmdMerger = installDirectory + "/bin/Donatello default_output_msa.fasta  outputMerger"
			subprocessLauncher(cmdPOA)
			subprocessLauncher(cmdMerger)
		if soft is not None:
			cmdMv = "mv outputMerger msa_" + soft + ".fa"
		else:
			cmdMv = "mv outputMerger msa.fa"
		subprocess.check_output(['bash','-c', cmdMv])
		cmdRM = "rm out*"
		subprocess.check_output(['bash','-c', cmdRM])
