# mlml - estimation of hydroxymethylation and methylation levels

## Synopsis
```
$ dnmtools mlml [OPTIONS] -u <input-bs.meth> -m <input-oxbs.meth> -h <input-tab.meth>
```

## Description

If you are interested in estimating hydroxymethylation level and have
any two of Tet-Assisted Bisulfite sequencing (TAB-seq), oxidative
bisulfite sequencing (oxBS-seq) and BS-seq data available, you can use
`mlml`  to perform consistent and simultaneous estimation.

The input file format could be the default
[counts](../counts) output format,
or BED format file with 6 columns as the example below:

```txt
chr1  3001345  3001346  CpG:9  0.777777777778  +
```

Here the fourth column indicates that this site is a CpG site, and the
number of reads covering this site is 9. The fifth column is the
methylation level of the CpG site, ranging from 0 to 1. Note that all
input files must be sorted. Assume you have three input files ready:
`input-bs-seq.meth`, `input-oxbs-seq.meth` and `input-tab-seq.meth`.
The following command will take all the inputs:

```shell
$ dnmtools mlml -v -u input-bs-seq.meth -m input-oxbs-seq.meth -h input-tab-seq.meth -o result.txt
```

If only two types of input are available, e.g. `input-bs-seq.meth and
`input-oxbs-seq.meth`, then use the following command:

```shell
$ dnmtools mlml -u input-bs-seq.meth -m input-oxbs-seq-seq.meth -o result.txt
```

In some cases, you might want to specify the convergence tolerance for
EM algorithm. This can be done through `-t` option.  For example:

```shell
$ dnmtools mlml -u input-bs-seq.meth -m input-oxbs-seq.meth -o result.txt -t 1e-2
```

This command will make the iteration process stop when the difference
of estimation between two iterations is less than 10^−2 . Thea value
format can be scientific notation, e.g. 1e-5, or float number, e.g.
0.00001.

The output of `mlml` is tab-delimited format. Here is an example:

```txt
chr11   15      16      0.166667        0.19697 0.636364        0
chr12   11      12      0.222222        0       0.777778        2
```

The columns are chromosome name, start position, end position, 5-mC
level, 5-hmC level, unmethylated level and number of conflicts. To
calculate the last column, a binomial test is performed for each input
methylation level (can be 2 or 3 in total depending on parameters). If
the estimated methylation level falls out of the confidence interval
calculated from input coverage and methylation level, then such event
is counted as one conflict. It is recommended to filter estimation
results based on the number of conflicts; if more conflicts happens on
one site then it is possible that information from such site is not
reliable.

## Options

```txt
 -o, -output
```
output file (default: STDOUT)

```txt
 -u, -bsseq
```
input BS-seq methcounts file
```txt
 -h, -tabseq
```
input TAB-seq methcounts file
```txt
 -m, -oxbsseq
```
input oxBS-seq methcounts file
```txt
 -t, -tolerance
```
EM convergence threshold (default: 0.000000)
```txt
 -a, -alpha
```
significance level of binomial test for each site (default: 0.050000)
```txt
 -H, -outh
```
hmC pseudo methcount output file (default: null)
```txt
 -M, -outm
```
mC pseudo methcount output file (default: null)
```txt
 -v, -verbose
```
print more run info to STDERR while the program is running.



