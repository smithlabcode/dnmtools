# merge - Combine counts files

## Synopsis
```
$ dnmtools merge [OPTIONS] <file1.meth> <file2.meth> ...
```

## Description

When working with a BS-seq project with multiple replicates, you may
first produce a methcounts output file for each replicate individually
and assess the reproducibility of the methylation result by comparing
different replicates. The `merge` program is used to merge
the those individual methcounts file to produce a single estimate that
has higher coverage. Suppose you have the three methcounts files from
three different biological replicates: `R1.meth`, `R2.meth` and
`R3.meth`.  To merge those individual methcounts files, execute

```
$ dnmtools merge R1.meth R2.meth R3.meth
```

The merge-methcounts program assumes that the input files end with an
empty line. The program can handle an arbitrary number of files, even
empty ones, that can have different numbers of lines.

## Options

```
 -o, -output
```
output file as [counts](../../analysis/counts) format (default: stdout)

```
 -h, -header
```
Print a header given by the input string at the top of the file
(ignored for tabular)

```
 -t, -tabular
```
Outputs the file as a table

```
 -f, -fractional
```
output fractions (requires `-tabular`)

```
 -r, -reads
```
minimum number of reads required when using the `-f` flag (default: 1)

```
 -v, -verbose
```
print more run info to STDERR while the program is running.


