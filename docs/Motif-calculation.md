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

- sample.cov.csv: The breadth of coverage with a certain depth of reads (e.g., 10x; 10x genome coverage / genome size).
  | _MAG ID_ | MAG 1 | MAG 1 | MAG 1 | MAG 3 ...|
  |:-- | :------: | :------: | :------: | :------: |
  | ___Sample 1___ | 0.98 | 0.5 | 0.3 | 0.75 ...|

- coverage.total.csv: Aggregation of sample.cov data tables. 
  | _MAG ID_ | MAG 1 | MAG 1 | MAG 1 | MAG 3 ...|
  |:-- | :------: | :------: | :------: | :------: |
  | ___Sample 1___ | 0.98 | 0.5 | 0.3 | 0.75 ...|
  | ___Sample 2___ | 0.92 | 0.4 | 0.82 | 0.55 ...|

- sample.motif.csv: The number of motif in the genome region with a certain depth.
  | _Motifs_ | GGAGG | GGCGG | GGAGA | GGCGA ...|
  |:-- | :------: | :------: | :------: | :------: |
  | ___MAG 1___ | 32 | 12 | 10 | 23 ...|
  | ___MAG 2___ | 20 | 11 | 12 | 34 ...|
  
- sample.methyl.csv: The number of methylated motif in the genome region with a certain depth.
  | _Motifs_ | GGAGG | GGCGG | GGAGA | GGCGA ...|
  |:-- | :------: | :------: | :------: | :------: |
  | ___MAG 1___ | 13 | 10 | 8 | 3 ...|
  | ___MAG 2___ | 10 | 9 | 10 | 33 ...|
  
- sample.frac.csv: The fraction table between _sample.methyl_ and _sample.motif_ (fraction = 100 * # methylated/ # motif).
  | _Motifs_ | GGAGG | GGCGG | GGAGA | GGCGA ...|
  |:-- | :------: | :------: | :------: | :------: |
  | ___MAG 1___ | 40.63 | 83.33 | 80.00 | 13.04 ...|
  | ___MAG 2___ | 50.00 | 81.81 | 83.33 | 97.06 ...|

- methylome.total.csv: Aggregation of sample.frac data tables. 
  | _MAG ID_ | Samples | Coverage | GGAGG | GGCGG | GGAGA ...|
  |:-- | :------: | :------: | :------: | :------: | :------: |
  | ___MAG 1___ | ___Sample 1___ | 0.98 | 40.63 | 83.33 | 80.00 ...|
  | ___MAG 1___ | ___Sample 2___ | 0.62 | 10.14 | 80.28 | 72.11 ...|
