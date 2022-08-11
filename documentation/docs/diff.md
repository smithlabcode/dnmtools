# diff - Compute methylation difference probabilities

## Synopsis
```
$ dnmtools diff [OPTIONS] <input-a.meth> <input-b.meth>
```

## Description

Suppose that we want to compare two methylomes: `input-a.meth` and
`input-b.meth`. We start out by calculating the differential
methylation score of each CpG site using the methdiff program:

```
$ dnmtools diff -o output.methdiff input-a.meth input-b.meth
```

Here are the first few lines of the output:

```
chr1    3000826 +       CpG     0.609908        16       7      21      11
chr1    3001006 +       CpG     0.874119        21      18      15      22
chr1    3001017 +       CpG     0.888384        20      19      15      25
chr1    3001276 +       CpG     0.010825         3      20      12      16
```

The first four columns are the same as the methcounts input. The fifth
column gives the probability that the methylation level at each given
site is lower in `input-a.meth` than `input-b.meth`. (For the other
direction, you can either swap the order of the two input files or
just subtract the probability from 1.0.) The method used to calculate
this probability is detailed in Altham (1971) [1], and can be thought
of as a one-directional version of Fisher's exact test. The remaining
columns in the output give the number of methylated reads of each CpG
in `input-a.meth`, number of unmethylated reads in `input-a.meth`,
number of methylated reads in `input-b.meth`, and number of
unmethylated reads in `input-b.meth`, respectively.

The methdiff output is used as input for the [dmr](../dmr) program.

## Options

```
 -p, -pseudo
```
Add a pseudocount to inputs (default: 1)

```
 -A, -nonzero-only
```
process only sites with coveage in both samples

```
 -o, -out
```
output file (default: STDOUT)

```
 -v, -verbose
```
print more run info to STDERR while the program is running.

