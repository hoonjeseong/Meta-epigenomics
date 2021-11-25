## __Meta-epigenomic analysis__

### __Table of contents__

- [Overview](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/README.md#overview)

- [Prerequisites](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/README.md#prerequisites)

### __Overview__

This toolkit helps detect DNA methylation patterns across metagenomic samples.

![0  Meta-Epi-min (1)](https://user-images.githubusercontent.com/39515472/143376782-f68a5aff-681a-4fc2-ad26-4c01439521b9.png)

### __Prerequisites__

- The metagenome-assembled genomes (MAGs) must be prepared prior to analysis as this is based on genome-centric metagenomics.
- This analysis requires the results of the base modification (ipdSummary) results from PacBio sequencing (https://libraries.io/github/ben-lerch/BaseMod-3.0). 
&larr; [Commands what I used for the ipdSummary is available here.](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/docs/tutorial-ipdSummary.md)

- This analysis requires _python 3.6 >=_ and the absolute paths of _bedtools_ and _samtools_ must be written in the _program.txt_ file as follows:

```
bedtools:[/usr/bin/bedtools] #(tested by bedtools v2.25.0)
samtools:[/usr/bin/samtools] #(tested by htslib 1.8)
 ```
 
- Required python library

Using the command of _pip3 install [somthing]_

`itertools, optparse, shutil, pathlib, gzip, pickle, tqdm, biopython`

Then, please copy this git

`git clone https://github.com/hoonjeseong/Meta-epigenome_analysis.git`

#The purpose of this snakemake workflow is to obtain high-quality metagenome-assembled genomes (MAGs) from previously generated assemblies. 


