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

#todo
def launchRscripts():
	pass


def generateLatexFigures( outDir, outputPDFName, filesDict):
	content = r'''\documentclass{article}
	\usepackage{graphicx}
	\begin{document}
	\section{Recall, precision of the correction and correct base rate in reads}
	\end{document} '''
	print(outDir + "/" + outputPDFName +'.tex')
	with open(outDir + "/" + outputPDFName +'.tex','w') as f:
		f.write(content%filesDict)
	proc = subprocess.Popen(['pdflatex', '-output-directory', outDir, outputPDFName + ".tex"])
	proc.communicate()


def generateResults(outDir):
	filesDict = {"recall": "recall.png"}
	generateLatexFigures(outDir, "summary", filesDict)
