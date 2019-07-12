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



def launchRscripts(installDirectory, soft, outDir):
	# recall and precision figure
	if soft is not None:
		if checkIfFile( outDir + "/" + soft + "_per_read_metrics.txt"):
			cmdRecallPrecision = "Rscript " + installDirectory + "/Rscripts/plot_recall_precision_correctrate.R " + outDir + "/" + soft + "_per_read_metrics.txt " + outDir
			subprocessLauncher(cmdRecallPrecision)
	else:
		if checkIfFile( outDir + "/per_read_metrics.txt"):
			cmdRecallPrecision = "Rscript " + installDirectory + "/Rscripts/plot_recall_precision_correctrate.R " + outDir + "/per_read_metrics.txt " + outDir
			subprocessLauncher(cmdRecallPrecision)

	# sizes distribution
	if soft is not None:
		if checkIfFile( outDir + "/" + soft + "_read_size_distribution.txt"):
			cmdSizesDistr = "Rscript " + installDirectory + "/Rscripts/plot_distribution_sizes.R " + outDir + "/" + soft + "_read_size_distribution.txt " + outDir
			subprocessLauncher(cmdSizesDistr)
	else:
		if checkIfFile( outDir + "/read_size_distribution.txt"):
			cmdSizesDistr = "Rscript " + installDirectory + "/Rscripts/plot_distribution_sizes.R " + outDir + "/read_size_distribution.txt " + outDir
			subprocessLauncher(cmdSizesDistr)



def generateLatexFigures( outDir, outputPDFName, filesDict, remap, assemble ):
	content = r'''\documentclass{article}
	\usepackage{graphicx}
	\usepackage{multirow}

	\begin{document}

	\section{Summary}
	\begin{figure*}[ht!]
	\begin{tabular}{|l|c|}
	\hline
	Assessed reads & %(nbReads)s \\ \hline
	Throughput & %(throughput)s \\ \hline\hline

	Recall (computed on corrected bases)& %(meanRecall)s  \\ \hline
	Precision (computed on corrected bases)& %(meanPrecision)s  \\ \hline
	Average correct bases rate  (computed on whole read)& %(meanCorrectBaseRate)s \\ \hline
	Overall error rate (computed on whole read)& %(errorRate)s \\ \hline\hline

	Number of trimmed/split reads & %(numberReadSplit)s \\ \hline
    Mean missing size in trimmed/split reads & %(meanMissingSize)s \\ \hline
    Number of reads over-corrected by extension & %(numberReadExtended)s \\ \hline
    Mean extension size in over-correcte reads & %(meanExtensionSize)s \\ \hline


    Number of corrected reads which length & \multirow{2}{*}{%(smallReads)s} \\
  is $<$ %(minLength)s \%% of the original read & \\ \hline\hline

  	Number of very low quality corrected reads & %(wronglyCorReads)s \\ \hline\hline

	Ratio of homopolymer sizes in corrected vs reference &  %(homoRatio)s \\ \hline\hline

	\%% GC in reference reads & %(GCRef)s \\ \hline
    \%% GC in corrected reads & %(GCCorr)s \\ \hline

	\end{tabular}
	\end{figure*}

	\begin{figure*}[ht!]
	\begin{tabular}{|l|c|c|}
	\hline
	 &Uncorrected & Corrected  \\ \hline
	Insertions & %(insU)s& %(insC)s  \\ \hline
	Deletions &  %(delU)s& %(delC)s  \\ \hline
	Substitutions &  %(subsU)s& %(subsC)s  \\ \hline
	\end{tabular}
	\end{figure*}
	'''

	if remap:
		content += r'''
	    \section{Corrected reads remapping on genome}
		\begin{figure*}[ht!]
		\begin{tabular}{|l|c|}
		\hline
		Average identity (\%%) & %(averageId)s  \\ \hline
		Genome coverage (\%%) & %(genomeCov)s  \\ \hline
		\end{tabular}
		\end{figure*}
		'''

	if assemble:
		content += r'''
		\section{Corrected reads assembly metrics}
		\begin{figure*}[ht!]
		\begin{tabular}{|l|c|}
		\hline
		Contigs number & %(nbContigs)s  \\ \hline
		Aligned contigs number & %(nbAlContig)s  \\ \hline
		Breakpoints number & %(nbBreakpoints)s  \\ \hline
		NGA50 & %(NGA50)s  \\ \hline
		NGA75 & %(NGA75)s  \\ \hline
		\end{tabular}
		\end{figure*}
		'''

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
		\caption{\textbf{Size distributions after correction.} "Sequences" relate to each fasta sequence of the corrected file. "Reads" relate to corrected reads. In case of split reads, a "read" can be composed of two "sequences" for instance. This is why we report two different distributions. In case no read is split, we only report read length distribution.}
		\label{fig:size_distr}
		\end{figure}'''
	content += r'''
	\end{document} '''
	with open(outDir + "/" + outputPDFName +'.tex','w') as f:
		f.write(content%filesDict)
	proc = subprocess.Popen(['pdflatex', '-output-directory', outDir, outputPDFName + ".tex"], stdout = DEVNULL, stderr = DEVNULL).communicate()
	#~ proc.communicate()



def generateResults(outDir, installDirectory, soft, nbReads, throughput, recall, precision, correctBaseRate, errorRate, numberSplit, meanMissing, numberExtended, meanExtension, percentGCRef, percentGCCorr, smallReads, wronglyCorReads, minLength, indelsubsUncorr, indelsubsCorr, avId, cov, nbContigs, nbAlContig, nbBreakpoints, NGA50, NGA75,  remap, assemble, homoRatio ):
	filesDict = {"recall_precision": outDir + "/plot_recall_precision.png", "size_distribution": outDir + "/plot_size_distribution.png", "nbReads": nbReads, "throughput": throughput, "meanPrecision": precision, "meanRecall": recall, "meanCorrectBaseRate": correctBaseRate, "errorRate": errorRate, "numberReadSplit": numberSplit, "meanMissingSize": meanMissing, "numberReadExtended": numberExtended, "meanExtensionSize": meanExtension, "GCRef": str(percentGCRef), "GCCorr": str(percentGCCorr), "smallReads": smallReads, "wronglyCorReads": wronglyCorReads, "minLength": minLength, "insC": indelsubsCorr[0], "delC": indelsubsCorr[1], "subsC": indelsubsCorr[2], "insU": indelsubsUncorr[0],"delU": indelsubsUncorr[1], "subsU": indelsubsUncorr[2], "averageId" : avId, "genomeCov": cov, "nbContigs": nbContigs, "nbAlContig" : nbAlContig, "nbBreakpoints": nbBreakpoints, "NGA50": NGA50, "NGA75": NGA75, "homoRatio": homoRatio}
	launchRscripts(installDirectory, soft, outDir)
	generateLatexFigures(outDir, "summary", filesDict, remap, assemble)



