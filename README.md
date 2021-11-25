## __The overview of the meta-epigenomic analysis__

#### - Overview

This toolkit helps detect DNA methylation patterns across metagenomic samples.

#### - Table of contents

#### - Prerequisites

This analysis requires the results of the base modification (ipdSummary) results from PacBio sequencing (https://libraries.io/github/ben-lerch/BaseMod-3.0). Please check the SMRTLink manual for recent chemicals and sequencers (e.g. Sequel II). I used pbbioconda (https://github.com/PacificBiosciences/pbbioconda) command to align raw reads, merge the multiple sequence result from cells, and detect base modifications.

Where reference.fasta contains your reference sequences.

Where subreads.bam contains your unaligned reads.

Where aligned_subreads.bam is where your aligned reads will be stored.

#The purpose of this snakemake workflow is to obtain high-quality metagenome-assembled genomes (MAGs) from previously generated assemblies. 

![0  Meta-Epi-min](https://user-images.githubusercontent.com/39515472/143327436-b50b0f01-ba70-4736-b5ba-e4e37bb1934e.png)



