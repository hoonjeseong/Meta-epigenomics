## __The overview of the meta-epigenomic analysis__

#### - Overview

This toolkit helps detect DNA methylation patterns across metagenomic samples.

#### - Table of contents

#### - Prerequisites

This analysis requires the results of the base modification (ipdSummary) results from PacBio sequencing (https://libraries.io/github/ben-lerch/BaseMod-3.0). Please check the SMRTLink manual for recent chemical and sequencer (e.g. Sequel II). I used pbbioconda command for alignment of raw reads, merging the multiple sequence result from cells, and detection of base modifications (https://github.com/PacificBiosciences/pbbioconda).

Where reference.fasta contains your reference sequences.

Where subreads.bam contains your unaligned reads.

Where aligned_subreads.bam is where your aligned reads will be stored.

#The purpose of this snakemake workflow is to obtain high-quality metagenome-assembled genomes (MAGs) from previously generated assemblies. 

![Meta-Epi](https://user-images.githubusercontent.com/39515472/143149711-66a11b4e-cdfa-45fe-b541-3cf3f8788fbc.png)


