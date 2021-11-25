## __Meta-epigenomic analysis__

### __Table of contents__

- [Overview](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/README.md#overview)
- [Prerequisites](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/README.md#prerequisites)
- [Getting MTase and motif information from REBASE](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/utils/MTase_REBASE.md)

### __Overview__

This toolkit helps detect DNA methylation patterns across metagenomic samples.

![0  Meta-Epi-min (1)](https://user-images.githubusercontent.com/39515472/143376782-f68a5aff-681a-4fc2-ad26-4c01439521b9.png)

### __Prerequisites__

- The metagenome-assembled genomes (MAGs) must be prepared prior to analysis as this is based on genome-centric metagenomics.
- This analysis requires the results of the base modification (ipdSummary) results from PacBio sequencing (https://libraries.io/github/ben-lerch/BaseMod-3.0). 
&larr; [Commands what I used for the ipdSummary is available here.](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/docs/tutorial-ipdSummary.md)

- The absolute paths of _bedtools_ and _samtools_ must be written in the _program.txt_ file as follows:

```
bedtools:[/usr/bin/bedtools] #(tested by bedtools v2.25.0)
samtools:[/usr/bin/samtools] #(tested by htslib 1.8)
 ```
 
- Required _python 3.6 >=_ and libraries

Using the command of _pip3 install [somthing]_

`itertools, optparse, shutil, pathlib, gzip, pickle, tqdm, biopython`

Then, please copy this git

`git clone https://github.com/hoonjeseong/Meta-epigenome_analysis.git`

- MTases and their motif information

If you want to get the motif of DNA methylation patterns reported in REBASE, then I suggest you access the analysis in [Getting MTase and motif information from REBASE](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/utils/MTase_REBASE.md). It provides motif information recognized by methyltransferases (MTases) in the [REBASE database](http://rebase.neb.com/rebase/rebase.html), for example in the form of G_A_NTC/G_A_NTC (forward/reverse strand).

### __Simple example__

#The purpose of this snakemake workflow is to obtain high-quality metagenome-assembled genomes (MAGs) from previously generated assemblies. 


