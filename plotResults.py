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
from utils import *



try:
    from subprocess import DEVNULL # py3k
except ImportError:
    import os
    DEVNULL = open(os.devnull, 'wb')


#todo
def launchRscripts(installDirectory, soft):
	# recall and precision figure
	#~ print("uuu", installDirectory)
	if soft is not None:
		if checkIfFile( installDirectory + "/" + soft + "_per_read_metrics.txt"):
			cmdRecallPrecision = "Rscript " + installDirectory + "/Rscripts/plot_recall_precision_correctrate.R " + installDirectory + "/" + soft + "_per_read_metrics.txt ."
			subprocessLauncher(cmdRecallPrecision)
	else:
		if checkIfFile( installDirectory + "/per_read_metrics.txt"):
			cmdRecallPrecision = "Rscript " + installDirectory + "/Rscripts/plot_recall_precision_correctrate.R " + installDirectory + "/per_read_metrics.txt ."
			subprocessLauncher(cmdRecallPrecision)
	
	# sizes distribution
	if soft is not None:
		if checkIfFile( installDirectory + "/" + soft + "read_size_distribution.txt"):
			cmdSizesDistr = "Rscript " + installDirectory + "/Rscripts/plot_distribution_sizes.R " + installDirectory + "/" + soft + "_read_size_distribution.txt ."
			subprocessLauncher(cmdSizesDistr)
	else:
		if checkIfFile( installDirectory + "/read_size_distribution.txt"):
			cmdSizesDistr = "Rscript " + installDirectory + "/Rscripts/plot_distribution_sizes.R " + installDirectory + "/read_size_distribution.txt ."
			subprocessLauncher(cmdSizesDistr)

def generateLatexFigures( outDir, outputPDFName, filesDict):
	content = r'''\documentclass{article}
	\usepackage{graphicx}
	\usepackage{multirow}

	\begin{document}

	\section{Summary}
	\begin{figure*}[ht!]
	\begin{tabular}{|l|c|} 
	\hline
	Mean recall & %(meanRecall)s  \\ \hline
	Mean precision & %(meanPrecision)s  \\ \hline
	Mean correct bases rate & %(meanCorrectBaseRate)s \\ \hline
   Number of trimmed/split reads & %(numberReadSplit)s \\ \hline
      Mean missing size in trimmed/split reads & %(meanMissingSize)s \\ \hline
    \%% GC in reference reads & %(GCRef)s \\ \hline
    \%% GC in corrected reads & %(GCCorr)s \\ \hline
    Number of corrected reads which length & \multirow{2}{*}{%(smallReads)s} \\
  is $<$ 10.0 \%% of the original read & \\ \hline
	\end{tabular}
	\end{figure*}


	\begin{figure*}[ht!]
	\begin{tabular}{|l|c|c|} 
	\hline
	 &uncorrected & corrected  \\ \hline
	Insertions & %(insU)s& %(insC)s  \\ \hline
	Deletions &  %(delU)s& %(delC)s  \\ \hline
	Substitutions &  %(subsU)s& %(subsC)s  \\ \hline
	\end{tabular}
	\end{figure*}'''

	
	#filesDict = {"recall_precision": installDirectory + "/plot_recall_precision.png", "size_distribution": installDirectory + "/plot_size_distribution.png", "meanRecall": recall, "meanPrecision": precision, "meanCorrectBaseRate": correctBaseRate, "numberReadSplit": numberSplit, "meanMissingSize": meanMissing, "GCRef": percentGCRef, "GCCorr": percentGCCorr, "smallReads": smallReads}

	#add plot recall precision
	if checkIfFile(filesDict["recall_precision"]):
		content += r'''
		\section{Recall, precision of the correction and correct base rate in reads}
		\begin{figure}[ht!]
		\centering\includegraphics[width=0.7\textwidth]{%(recall_precision)s}
		\caption{\textbf{Recall and precision computed for each read after correction.} Recall is the rate of correct bases out of bases that needed correction in the original read, precision is the rate of correct bases out of the number of bases that were modified by the correction method.}
		\label{fig:recall_prec}
		\end{figure}'''

	# add length distribution
	if checkIfFile(filesDict["size_distribution"]):
		content += r'''
		\begin{figure}[ht!]
		\centering\includegraphics[width=0.7\textwidth]{%(size_distribution)s}
		\caption{\textbf{Size distribution of reads before and after correction.}}
		\label{fig:size_distr}
		\end{figure}'''
	content += r'''\end{document} '''
	with open(outDir + "/" + outputPDFName +'.tex','w') as f:
		f.write(content%filesDict)
	proc = subprocess.Popen(['pdflatex', '-output-directory', outDir, outputPDFName + ".tex"], stdout = DEVNULL, stderr = DEVNULL).communicate()
	#~ proc.communicate()


def generateResults(outDir, installDirectory, soft, recall, precision, correctBaseRate, numberSplit, meanMissing, percentGCRef, percentGCCorr, smallReads, indelsubsUncorr, indelsubsCorr, avId, cov, nbContigs, nbAlContig, nbBreakpoints, NG50, NG75 ):
	filesDict = {"recall_precision": installDirectory + "/plot_recall_precision.png", "size_distribution": installDirectory + "/plot_size_distribution.png", "meanRecall": recall, "meanPrecision": precision, "meanCorrectBaseRate": correctBaseRate, "numberReadSplit": numberSplit, "meanMissingSize": meanMissing, "GCRef": str(percentGCRef), "GCCorr": str(percentGCCorr), "smallReads": smallReads, "insC": indelsubsCorr[0], "delC": indelsubsCorr[1], "subsC": indelsubsCorr[2], "insU": indelsubsUncorr[0],"delU": indelsubsUncorr[1], "subsU": indelsubsUncorr[2]}
	launchRscripts(installDirectory, soft)
	generateLatexFigures(outDir, "summary", filesDict)


#Recall  0.9997 
#Precision  1.0 
#Correct bases rate  0.99998 
#Number of trimmed/split reads  0 
#Mean missing size in trimmed/split reads  0
#%GC in reference reads: 36.1 %GC in corrected reads: 18.5
#Number of corrected reads which length is < 10.0 % of the original read: 0
