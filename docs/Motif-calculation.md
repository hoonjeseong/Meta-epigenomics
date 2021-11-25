## __Motif calculation__

![2 Meta-Epi_p (1)](https://user-images.githubusercontent.com/39515472/143390947-9b998ed5-e063-4f8a-877b-53ca5cf09266.png)

_The purpose of this script is to obtain methylated frequency of each motif on MAGs across metagenomic samples._
- Search motifs and methylated motifs from REBASE (___mtase_pickle___ file from the [document](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/8bd14c73bd7b95333fe8fa2be10bf505c827ed57/utils/MTase_REBASE.md))  or _de novo_ ___5bp motifs___ (NN_A_NN or NN_C_NN; _A_ and _C_ means modified base of 6mA and 4mC which could be detected by PacBio system) on MAGs by each metagenomic sample.
- For simple normalization, only these motifs are counted in 10 dp mapped bredth regions, regardless of the strands in the MAG.
