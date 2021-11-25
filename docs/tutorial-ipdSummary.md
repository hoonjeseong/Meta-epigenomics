## __PacBio DNA modification detection methods__

### __Prerequisites__

  This analysis requires the results of the base modification (ipdSummary) results from PacBio sequencing (https://libraries.io/github/ben-lerch/BaseMod-3.0). I used pbbioconda (https://github.com/PacificBiosciences/pbbioconda) command to ___align___ raw reads, ___merge___ the multiple sequence result from cells, and ___detect___ base modifications.
  
  However, for sequence data using recent chemistries and sequencers (e.g. Sequel II) please check the SMRTLink manual before used this commands.

### __Command line__

- __Load the environment__

`source activate pbcore`


- __Alignment subreads to MAGs__

```
pbalign --concordant --hitPolicy=randombest --minAccuracy 70 --minLength 50 \
    --algorithmOptions="--minMatch 12 --bestn 10 --minPctIdentity 70.0" --nproc [#threads] \
    [input subreads.bam] [reference.fa] [aligned subreads.bam]
```
___#threads___ is number of threads you want to use.

___input subreads.bam___ is sequencing reads of input.

___reference.fa___ is reference genome fasta file which want to know (e.g., MAGs)

___aligned subreads.bam___ is the aligned result to `reference.fa`


- __Merge subreads from same sample__

If you have multiple subreads (PacBio cells) from a same sample for increasing sequecning throuput, merge the alignment result before calculate base modifications.

`pbmerge -o [merged sample.bam] [aligned *subreads.bam]`

___aligned *subreads.bam___ is aligned bam files you want to merge. Using wildcard is also allowd.

___merged sample.bam___ is a merged bam file from same sample.


- __Merge subreads from same sample__

```
ipdSummary [merged sample.bam] --reference [reference.fa] --gff [sample.gff] --csv [sample.csv] \
  --pvalue 0.001 --numWorkers [#threads] --identify m4C,m6A --minCoverage 3 --methylMinCov 10 \
  --methylFraction
```
___sample.gff___ is a DNA modification result.

___sample.csv___ is a detailed DNA modification result.

_sample.gff and sample.csv files are needed for further Meta-Epigenomic analysis_
