# multistat

## Synopsis
```
$ dnmtools multistat [OPTIONS] <intervals.bed> <probes.bed> <input.meth>
```

## Description

program to summarize methylation values according to genomic intervals
specified in a BED format file. The column headings must be the names
associated with intervals in a separate BED file. The methylation
table must have columns corresponding to sites in the genome.

## Options
```
 -o, -outfile
```
output file
```
 -v, -verbose
```
print more run info
```
 -progress
```
report progress
```
 -name-by-interval
```
name features by interval in output
```
 -e, -empty
```
report empty intervals
