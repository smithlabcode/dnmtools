# entropy - Computing methylation entropy

## Synopsis
```shell
$ dnmtools entropy [OPTIONS] <genome.fa> <input.epiread>
```
## Description
The concept of methylation entropy was introduced into epigenetics
study to characterize the randomness of methylation patterns over
several consecutive CpG sites (Xie et al, 2011). The `methentropy`
program processes epireads and calculates the methylation entropy
value in sliding windows of specified number of CpGs. Two input files
are required.

 * (1) either a genome in FASTA format or a directory containing FASTA
   chromosome files files

 * (2) an epiread file as produced by
   [states](../states) program. The input epiread file
   needs to be sorted, first by chromosome, then by position. It can
    be done with the following command.

```shell
$ LC_ALL=C sort -k1,1 -k2,2g input.epiread -o input-sorted.epiread
```

Use the `-w` option to specify the desired number of CpGs in the
sliding window; if unspecified, the default value is 4. In cases where
symmetric patterns are considered the same, specify option -F, this
will cause the majority state in each epiread to be forced into
"methylated", and the minority to "unmethylated". The processed
epireads will then be used for entropy calculation. To run the
program, type command:
```shell
$ dnmtools entropy -w 5 -v -o output.meth /path/to/genome.fa input-sorted.epiread
```

The output format is the same as [counts](../counts)
output. The first 3 columns indicate the genomic location of the
center CpG in each sliding window, the 5th column contains the entropy
values, and the 6th column shows the number of reads used for each
sliding window.  Below is an output example.

```txt
chr1    483     +       CpG     2.33914 27
chr1    488     +       CpG     2.05298 23
chr1    492     +       CpG     1.4622  24
chr1    496     +       CpG     1.8784  35
```

## Options
```txt
 -w, -window
```
number of CpGs in sliding window (default: 4)
```txt
 -F, -flip
```
flip read majority state to meth
```txt
 -o, -output
```
Name of output file (default: STDOUT)
```txt
 -v, -verbose
```
print more run info to STDERR while the program is running.

