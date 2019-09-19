Reproduce results of ELECTOR's paper
===================================

First, you need to download, install, and have the following tools (and their dependencies) available in your $PATH.

1) For read simulation :

* ART (https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
* SimLoRD (https://bitbucket.org/genomeinformatics/simlord/src/master/)

2) For error correction (you do not need all the tools, only the ones you wish to use):

* Canu (https://github.com/marbl/canu)
* Daccord (https://github.com/gt1/daccord)
* HG-CoLoR (https://github.com/morispi/HG-CoLoR)
* LoRDEC (http://www.atgc-montpellier.fr/lordec/)
* MECAT (https://github.com/xiaochuanle/MECAT)

Scripts are provided to recreate the experiments on the following organisms:

* E. coli
* S. cerevisiae
* C. elegans

You first need to simulate the short reads, using the scripts in the ``shortReads`` subdirectory.
Then, simulate the long reads, using the scripts in the ``longReads`` subdirectory.
Pick the correction method you wish to use, and launch it using the scripts in the ``correction/correctionMethod`` directory.
Finally, evaluate the quality of the correction, using the scripts in the ``evaluation/correctionMethod`` directory.

For instance, if you wish to run the evaluation of MECAT on the E. coli dataset, proceed as follows:

	cd shortReads/
	./simEcoli.sh
	cd ../
	cd longReads/
	./simEcoli.sh
	cd ../
	cd correction/
	cd MECAT/
	./corEcoli.sh
	cd ../../
	cd evaluation/
	cd MECAT/
	./evalEcoli.sh
