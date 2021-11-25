## __The overview of the meta-epigenomic analysis__

#### - Overview

This toolkit helps detect DNA methylation patterns across metagenomic samples.

#### - Table of contents

#### - Prerequisites

This analysis requires the results of the base modification (ipdSummary) results from PacBio sequencing (https://libraries.io/github/ben-lerch/BaseMod-3.0). 
&larr; [Commands what I used for the ipdSummary is available here.](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/docs/tutorial-ipdSummary.md)

This analysis requires _python 3.6 >=_ and the absolute paths of _bedtools_ and _samtools_ must be written in the _program.txt_ file as follows:

```
bedtools:[/usr/bin/bedtools] #(tested by bedtools v2.25.0)
samtools:[/usr/bin/samtools] #(tested by htslib 1.8)
 ```
 
_required python library_

using the command of _pip3 install [somthing]_

```
- itertools
- optparse
- shutil
- pathlib
- gzip
- pickle
- tqdm
- biopython
```

`git clone https://github.com/hoonjeseong/Meta-epigenome_analysis.git`

#The purpose of this snakemake workflow is to obtain high-quality metagenome-assembled genomes (MAGs) from previously generated assemblies. 

![0  Meta-Epi-min](https://user-images.githubusercontent.com/39515472/143327436-b50b0f01-ba70-4736-b5ba-e4e37bb1934e.png)



