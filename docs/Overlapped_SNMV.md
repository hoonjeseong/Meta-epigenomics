## __Overlapped ingle nucleotide methylated variation (SNMV)__

![3 Meta-Epi_p (1)](https://user-images.githubusercontent.com/39515472/143391267-e5f5a76a-6eaa-45db-8a36-a00fd18ae5df.png)

_The purpose of this script is to calculate the fraction of specific methylated motif by each nucleotide position on MAGs across metagenomic samples._
- Search motifs and methylated motifs from REBASE (___mtase_pickle___ file from the [document](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/8bd14c73bd7b95333fe8fa2be10bf505c827ed57/utils/MTase_REBASE.md))  or _de novo_ ___5bp motifs___ (NN_A_NN or NN_C_NN; _A_ and _C_ means modified base of 6mA and 4mC which could be detected by PacBio system) on MAGs by each metagenomic sample.
- For simple normalization, only these motifs are counted in 10 dp mapped bredth regions, regardless of the strands in the MAG.
