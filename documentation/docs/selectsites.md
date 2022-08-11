# selectsites - Subsetting cytosines in methcounts files

## Synopsis
```shell
 $ dnmtools selectsites [OPTIONS] <regions.bed> <input.meth>
```

In many cases, we may be interested in analyzing only a subset of
cytosines or CpGs in a sample. Some instances of these cases including
calculating average methylation levels in (1) annotated regions, such
as promoter regions or repeats or (2) regions defined by the data
itself, such as HMRs or PMDs.

A possible solution to subset these regions is to convert the
methcounts file to BED format, intersect it with a BED file of the
regions of interest (using
[bedtools](https://bedtools.readthedocs.io)), then convert it back to
methcounts. The program selectsites simplifies these operations. It
takes a [counts](../counts) file and a set of intervals
represented as a BED file and produces a subset of the entries in the
methcounts file included in the BED regions. We can select entries in
input.meth contained in any inverval in `regions.bed` using the
following command.

```shell
$ dnmtools selectsites -o output.meth regions.bed input.meth
```

## Options

```txt
 -o, -output
```
Name of output file (default: STDOUT)

```txt
 -p, -preload
```
preload sites (use for large target intervals)

```txt
 -v, -verbose
```
print more run info to STDERR while the program is running.
