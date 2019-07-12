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

#!/usr/bin/env python3
import sys
import os
import shlex
import subprocess
from subprocess import Popen, PIPE, STDOUT
import re
import copy
import argparse
import glob



installDirectory = os.path.dirname(os.path.realpath(__file__))+"/../bin/"
dataDirectory = os.path.dirname(os.path.realpath(__file__))+"/../src/poa-graph/"



try:
    from subprocess import DEVNULL # py3k
except ImportError:
    import os
    DEVNULL = open(os.devnull, 'wb')



######### utils for warnings, user messages, errors #########
#warnings



def printWarningMsg(msg):
	print("[Warning] " + msg)



# to return if an error makes the run impossible
def dieToFatalError (msg):
  print("[FATAL ERROR] " + msg)
  sys.exit(1);



# check if file exists and is not empty
def checkIfFile(pathToFile):
	if not(os.path.exists(pathToFile) and os.path.getsize(pathToFile) > 0):
		return False
	return True



def checkReadFiles(readfiles):
	if readfiles is None:
		return True
	allFilesAreOK = True
	#~ for file in readfiles:
	if not os.path.isfile(readfiles):
		print("[ERROR] File \""+readfiles+"\" does not exist.")
		allFilesAreOK = False
	if not allFilesAreOK:
		dieToFatalError("One or more read files do not exist.")



######### utils for subprocess #########
# launch subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
	args = shlex.split(cmd)
	#~ p = subprocess.call(args, stdin = argstdin, stdout = argstdout, stderr = argstderr)
	p = subprocess.Popen(args, stdin = argstdin, stdout = DEVNULL, stderr = DEVNULL).communicate()
	return p



######### utils for sequence files #########
# find files with a regex
def getFiles(pathToFiles, name): #for instance name can be "*.txt"
	os.chdir(pathToFiles)
	listFiles = []
	for files in glob.glob(name):
		listFiles.append(files)
	return listFiles



# return number of reads in a fasta
def getFileReadNumber(fileName):
	cmdGrep = """grep ">" -c """ + fileName
	val = subprocess.check_output(['bash','-c', cmdGrep])
	return int(val.decode('ascii'))



def getCorrectedSequence(fileName, header):
	cmdGrep = """grep -A 1 """ + header + " " + fileName + "| tail -1"
	val = subprocess.check_output(['bash','-c', cmdGrep])
	return str(val.decode('ascii'))
