# multistat

## Synopsis
```shell
$ dnmtools multistat [OPTIONS] <intervals.bed> <input-tabular.tsv>
```

## Description

The `multistat` program is similar to [roi](../roi), but instead of creating a
BED file with averge methylation levels from a single counts file, it takes as
an input the output of [merge](../merge) with tabular format (i.e. using the
`-t` flag to make a data frame). In other words, `multistat` takes a data frame
as input and produces a data frame as output.

The input of `multistat` starts with a line with `2n` column names, with each
column name appearing sequentially twice. The file is then followed by a set of
lines containing `2n+1` elements. Each sample contains two columns. The first
column is the number of reads that cover the CpG in the sample, and the second
column is the number of CpGs that are methylated among the reads.

Here is a visual example of a file called `input-tabular.tsv` with four samples
(D083a, D083b, D091a and D091b):

```txt
  D083a       D083a       D083b       D083b       D091a       D091a       D091b       D091b
chr1:10468:+:CpG        3       0       2       0       2       0       1       0
chr1:10470:+:CpG        6       0       3       0       4       0       3       0
chr1:10483:+:CpG        7       0       3       0       5       0       3       1
chr1:10488:+:CpG        7       0       3       0       5       0       3       0
chr1:10492:+:CpG        7       0       2       0       5       0       3       0
chr1:10496:+:CpG        6       0       4       0       5       0       4       0
chr1:10524:+:CpG        6       2       4       0       7       0       5       1
chr1:10541:+:CpG        4       0       4       0       7       2       5       0
chr1:10562:+:CpG        3       0       3       0       6       0       4       0
chr1:10570:+:CpG        2       0       3       0       6       0       4       0
chr1:10576:+:CpG        2       0       3       0       6       0       4       0
```

Note that the output of `merge` with `-t` flag may contain suffixes `_R` and
`_M` on the column names (e.g. `D083a_R` and `D083a_M` corresponding to the
"Reads" and "Methylated" columns). You can remove these to make the input proper
by running

```shell
 $ sed -i '1s/_[MR]//g' input-tabular.tsv
```

`multistat` also requires an input BED file representing the genomic
intervals of interest. The regions must be sorted by chromosome, position and
strand. If they are not, you can add the `-s` flag to sort the file prior to
running the program. Note that for very large BED files, this may take a long
time. Given an input file `regions.bed`, you can sort it in one of the two
following ways:

```shell
 $ LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k6,6 -o regions.bed regions.bed
```

```shell
 $ bedtools sort -i regions.bed
```

Finally, to create a file `data-frame.tsv` with methylation levels (which can be
[weighted, unweighted or fractional](../levels) methylation), run

```shell
 $ dnmtools multistat -o data-frame.tsv regions.bed input-tabular.tsv
```

## Options

```txt
 -o, -output
```

Name of output file (default: STDOUT)

```txt
 -N, -numeric
```

print numeric values only (not NAs), guaranteeing that the output
contains as many rows as there are regions in the BED input.

```txt
 -L, -preload
```

Load all CpG sites

```txt
 -s, -sort
```

sort data if needed


```txt
 -l, -level
```

the level to report as score column in bed format output (w, u or f),
corresponding to weighted, unweighted or fractional methylation (default: w)

```txt
 -M, -more-levels
```

report more methylation information

```txt
 -v, -verbose
```

print more run info to STDERR
