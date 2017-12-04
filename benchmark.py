#!/usr/bin/env python3

"""*****************************************************************************
 *   Authors: Camille Marchet Antoine Limasset Pierre Morisse Lolita Lecompte
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



import sys
import os
import shlex, subprocess
from subprocess import Popen, PIPE, STDOUT
import re



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



def main():
	beg = time.time()
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	# Manage command line arguments
	parser = argparse.ArgumentParser(description="Benchmark for quality assessment of long reads correctors.")

	# Define allowed options
	#~ parser.add_argument('MAPPINGOUTPUT', action="store", help="output file from mapper")
	#~ parser.add_argument('-a', action="store", help="annotation file (.gtf) , not mandatory", dest="annotationFile",
					  #~ default=None)
	#~ parser.add_argument('-o', action="store", dest="outputF", default="events.txt",  help="outputFileName")
	#~ parser.add_argument('--del_IR', help="threshold to decipher if an event is a deletion or an intron retention. See "
									   #~ "model for more information. Default 50",
					  #~ action="store", dest="threshold", default=50)
	#~ parser.add_argument('--oo',  help="name for not regular output files. For instance if 'foo.txt' is given, the 3 not "
									#~ "regular output files will be stored in the same directory than the regular output"
									#~ "file, with the names 'exact_repeats_foo.txt', 'inexact_repeats_foo.txt', "
									#~ "'not_mapped_events_foo.txt'", action="store", dest="otherOutputsF",
					  #~ default="")
	#~ parser.add_argument('--counts', help="must match with the value of KisSPlice --counts option, default : 0",
					  #~ action="store", dest="countOption", default=0, choices=(0, 1, 2), type=int)
	
	# get options for this run
	#~ options = parser.parse_args()
	end = time.time()
	print "Run ends in {0} seconds.".format(str(round(end-beg, 2)))


		
if __name__ == '__main__':
	main()
