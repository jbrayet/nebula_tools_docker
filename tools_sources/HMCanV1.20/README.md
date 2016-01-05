HMCan is a tool to detect histone modifications in cancer samples. 

Author: Haitham Ashoor
Contact: haitham.ashoor@kaust.edu.sa,Valentina.Boeva@curie.fr

Compilation:
In order to compile HMCan on your system use the provided make file with HMCan
In order to complie the HMCan, write the following command:
> make 

Running HMCan:
you can run HMCan using the following command:
./HMCan <TargetFile> <ControlFile> <configuration file> <Name>

TargetFile: Alignment file for target ChIP-seq data. 
ControlFile: Alignment file for control ChIP-seq data.
Name: run name. Output files will have the name as prefix. 


Configuration file:
In configuration file you can provide parameters for HMCan algorithm, 
the description of the parameters is described below.

default parameters are provided in HMCan.config.txt.

Parameters description:

1)format: HMCan accepts BAM,SAM and BED alignment formats. In order to process BAM file, 
samtools should be installed on the system.

2)GCIndex: this file contains GC content and mapability scores for large regions of DNA. The
file should be formated in cnp format. We provide "GC_profile_100KbWindow_Mapp76_hg19.cnp"  and "GC_profile_100KbWindow_Mapp76_hg19_ensemble.cnp" assample files. 
Please set largeBinLength to 10000 when using this file.

3)genomePath: a path for the sequences of the genome under study, each file corresponds to
one chromosome. One can download hg19 chromosomes from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/

4)minLength: minimum fragment length used in the ChIP-seq experiment 

5)medLength: median fragment length used in the ChIP-seq experiment

6)maxLength: maximum fragment length used in the ChIP-seq experiment

7)smallBinLength: bin size to be used to construct signal profiles

8)largeBinLength: bin size to be used to calculate copy number, 
please note that should be the same for GCIndex 

9)pvalueThreshold: p-value threshold of Poisson's single side exact test 

10)mergeDistance: max distance to merge single peaks into region

11)iterationThreshold: Score threshold to remove peaks in the iteration stage

12)finalThreshold : Score threshold to report peaks or regions 

13)maxIter: maximum iterations for HMCan algorithm

14)printWig: a boolean option enables the user to print density in WIG files or not.

15)printPosterior: a boolean option enables the user to print the bins posterior probabilities in a WIG format.

16) PorsteriorProb: threshold for posterior probability to consider bin a peak state

17) blackListFile: a file that contains regions to be ignored by HMCan. It should be  a three columns file chr,start,end. Example for Human blacklist file is in data folder


Output files:
HMCan produces Three output files for histone peaks, regions, and whole genome density. Peaks and regions files are in BED format,
the last column of the file corresponds to the peak or region summit. Density files are in wig format. 


Citation:
If you use HMCan please cite HMCan: a method for detecting chromatin modifications in cancer samples using ChIP-seq data 
Haitham Ashoor; Aurelie Herault; Aurelie Kamoun; Francois Radvanyi; Vladimir B. Bajic; Emmanuel Barillot; Valentina Boeva
Bioinformatics 2013; doi: 10.1093/bioinformatics/btt524

NOTE: starting from V1.16 chromosome names in the GC index file should match exactly names in alignment and provided genome file. We provide two variants of GC index file. 
