# selectsites - get subsets of cytosines from counts files

## Synopsis
```shell
$ dnmtools selectsites [OPTIONS] <regions.bed> <input.counts>
```

## Description

In many cases, we may be interested in analyzing only a subset of
cytosines or CpGs in a sample. Some instances of these cases including
calculating average methylation levels in (1) annotated regions, such
as promoter regions or repeats or (2) regions defined by the data
itself, such as HMRs or PMDs.

A possible solution to subset these regions is to convert the counts file to
BED format, intersect it with a BED file of the regions of interest (using
[bedtools](https://bedtools.readthedocs.io)), then convert it back to
counts. The program selectsites simplifies these operations. It takes a
[counts](../counts) format file and a set of intervals in a BED file and
produces a subset of the entries in the counts file included in the BED
regions. We can select entries in `input.counts` contained in any inverval in
`regions.bed` using the following command.

```shell
$ dnmtools selectsites -o output.counts regions.bed input.counts
```

## Options

```txt
 -o, -output
```
Name of output file (default: STDOUT)

```txt
 -p, -preload
```
Preload sites (use for large target intervals).

```txt
 -v, -verbose
```
Print more run info to STDERR while the program is running.

```txt
 -d, -disk
```
Process sites on disk (fast if target intervals are few).

```txt
 -S, -summary
```
Write summary to this file.

```txt
 -z, -zip
```
The output file will be in gzip compressed format.

```txt
 -relaxed
```
Allow additional columns in the input file.
