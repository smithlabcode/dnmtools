# diff - compute methylation difference probabilities

## Synopsis
```console
$ dnmtools diff [OPTIONS] <input-a.meth> <input-b.meth>
```

## Description

Suppose that we want to compare two methylomes: `input-a.meth` and
`input-b.meth`. Both these files would have been produced by the
[counts](../counts) command. We start by calculating the differential
methylation score (probability) for each CpG site using the `diff`
command:

```console
$ dnmtools diff -o output.diff input-a.meth input-b.meth
```

Here are the first few lines of the output:

```txt
chr1    3000826   +     CpG     0.609908        16       7      21      11
chr1    3001006   +     CpG     0.874119        21      18      15      22
chr1    3001017   +     CpG     0.888384        20      19      15      25
chr1    3001276   +     CpG     0.010825         3      20      12      16
```

The first four columns are the same as the counts input files. The 5th
column gives the probability that the methylation level at each given
site is lower in `input-a.meth` than `input-b.meth`. (For the other
direction, you can either swap the order of the two input files or
just subtract the probability from 1.0.) The method used to calculate
this probability is explained by Altham (see reference below), and is
simply a one-directional version of Fisher's exact test. The remaining
columns in the output give the number of methylated reads of each CpG
in `input-a.meth`, number of unmethylated reads in `input-a.meth`,
number of methylated reads in `input-b.meth`, and number of
unmethylated reads in `input-b.meth`, respectively.

The two input files must be have all sites within a chromosomes
consecutive, have the same chromosome order, and have sites sorted in
increasing order within each chromosome. The order of chromosomes does
not matter (e.g., chr10 may precede chr2, or chr2 may precede chr10).

**Warning** the order of the samples/methylomes given as input, the
"a" and "b", matters. It is probably a good idea to include this order
in the output file name, for example as `output_a_lt_b.diff`.

The output from the `diff` command is used as input for the
[dmr](../dmr) program, but may also form the basis of visualization if
you want to plot differential methylation probabilities, for example
along the genome in a genome browser.

Reference:
```txt
Patricia M. E. Altham (1969)
Exact bayesian analysis of a 2x2 contingency table, and Fisher's "exact" significance test
Journal of the Royal Statistical Society, Series B (Methodological)
31(2):261-269
```

## Options

```txt
-p, -pseudo
```
The pseudocount to use (default: 1).

```txt
-A, -nonzero-only
```
Process only sites with coveage in both samples.

```txt
-o, -out
```
The name of the output file. If no file name is provided, the output
will be written to standard output. Due to the size of this output, a
file name should be specified unless the output will be piped to
another command or program.

```txt
-v, -verbose
```
Print more information while the command is running.
