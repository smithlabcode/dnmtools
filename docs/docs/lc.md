# lc - Count number of lines in a big file

## Synopsis
```
$ dnmtools lc <input-human-readable-file>
```

## Description

When working with next-generation sequencing data, researchers often
handle very large files, such as FASTQ files containing raw reads and
\*.sam files containing mapped reads. `lc_approx` is an auxiliary tool
designed to approximate the number of lines in a very large file by
counting the number of lines in a small, randomly chosen chunk from
the big file and scaling the estimate by file size. For example, in
order to estimate the number of reads in a FASTQ file `input.fq`, run
```
$ dnmtools lc input.fq
```
It will return the approximate number of lines in this file and by
dividing the above number by 4, you get the approximate number of
reads in that file. The lc approx can be hundreds of times faster than
the unix tool` wc -l`.

## Options

```
 -v, -verbose
```
print more run info to STDERR while the program is running.
```
 -n, -samples
```
number of samples
```
 -z, -size
```
sample size (bytes)

