## __Motif calculation__

![meta-epi2](https://user-images.githubusercontent.com/39515472/143435477-dbab5ea5-3a4c-4371-98a3-66ecd30b8058.png)

_The purpose of this script is to obtain methylated frequency of each motif on MAGs across metagenomic samples._
- Search motifs and methylated motifs from REBASE (___mtase_pickle___ file from the [document](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/utils/MTase_REBASE.md))  or palindromic ___5bp motifs___ (NN_A_NN or NN_C_NN; _A_ and _C_ means modified base of 6mA and 4mC which could be detected by PacBio system) on MAGs by each metagenomic sample.
- For simple normalization, only these motifs are counted in 10 dp mapped bredth regions, regardless of the strands in the MAG.

### __- Usage__
```
python Motif_calculation.py -b [bam; folder of bam files; Extension: bam] -i [folder of ipdSummary files; Extension: gff] -g [folder of MAGs; Extension: fa or fna] -o [output]
```

___folder of bam files___: _sample.bam_ file from [aligned and merged subreads](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/utils/MTase_REBASE.md) is required.

___folder of ipdSummary files___: _sample.gff_ file from [ipdSummary](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/utils/MTase_REBASE.md) is required.

___folder of MAGs___: The folder where fasta file of MAGs are located. 

___output___: The output folder 
  
This will create the following files in the following output folder: 

- sample.cov.csv: A log file

- sample.motif.csv: The assignments from NMF unmanipulated useful for identifying multicopy genes.

- sample.frac.csv: As above but discretised NMF predictions.

- coverage.total.csv: Predictions from the Gibbs sampler selecting run with maximum log posterior.

- methylome.total.csv: 
  
  | _MAG ID_ | MAG 1 | MAG 1 | MAG 1 | MAG 3 ...|
  |:-- | :------: | :------: | :------: | :------: |
  | ___Contig ID:base position:strand___ | MAG1_1:5:+ | MAG1_1:20:+ | MAG1_1:23:- | MAG1_3:115:+ ...|
  | ___Sample 1___ | 1.0 (methylation fraction) | 0.8 | 0.9 | 0.75 ...|
  | ___Sample 2___ | 0.3 | 0.9 | 0.7 | 0.95 ...|
  | ___Sample 3___ | 0.3 | 0.9 | 0.7 | 0.95 ...|
