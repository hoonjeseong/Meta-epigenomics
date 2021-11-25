## __Overlapped single nucleotide methylated variation (SNMV)__

![meta-epi3](https://user-images.githubusercontent.com/39515472/143435335-04425562-294b-4c47-a949-e67d4add8852.png)

_The purpose of this script is to calculate the fraction of specific methylated motif by each nucleotide position on MAGs across metagenomic samples._

- Before comparison of methylation fraction, you should find DNA methylated motifs of each MAG.
- DNA methylated motifs could be inferred by [finding the highest DNA methylated motif](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/docs/Motif-calculation.md) and [comparing MTase similarity with its recognition motif](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/utils/MTase_REBASE.md).

```
Usage: python Overlapped_SNMV.py -i [folder of ipdSummary files; Extension: csv and gff] -g [folder of MAGs; Extension: fa or fna] -m motif [fwd/rev] -o [output]
```

___folder of ipdSummary files___: _sample.csv_ and _sample.gff_ files from [ipdSummary](https://github.com/hoonjeseong/Meta-epigenome_analysis/blob/main/utils/MTase_REBASE.md) are required.

___folder of MAGs___: The folder where fasta file of MAGs are located. 

___motif (fwd/rev)___: Forward and reverse motifs are required for input. Example: __'G_A_NTC/G_A_NTC'__ (A nucleotide surrounded by two underscores indicates a DNA modification)

___output___: The output table.
