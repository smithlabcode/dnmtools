# allelic - Single-site ASM scoring

## Synopsis
```shell
$ dnmtools allelic [OPTIONS] <input.epiread>
```

## Description

The program `allelicmeth`  calculates allele specific methylation scores
for each CpG site. Input files should be the epiread files (.epiread
suffix) produced using [states](../states). In the output file, each row
represents a CpG pair made by any CpG and its previous CpG, the first
three columns indicate the positions of the CpG site, the fourth
column is the name including the number of reads covering the CpG
pair, the fifth column is the score for ASM, and the last four columns
record the number of reads of four different methylation combinations
of the CpG pair: methylated methylated (mm), methylated unmethylated
(mu), unmethylated methylated (um), or unmethylated unmethylated (uu).
The following command will calculate allele specific methylation
scores using the allelicmeth component of dnmtools:

```shell
$ dnmtools allelic -c /path/to/genome.fa -o output.allelic input.epiread
```

## Options

```txt
 -o, -output
```
output file name (default: STDOUT)
```txt
 -c, -chrom
```
FASTA file or directory of chromosomes containing FASTA files [required]

```txt
 -v, -verbose
```
print more run info to STDERR while the program is running.
