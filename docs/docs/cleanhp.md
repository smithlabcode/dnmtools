# cleanhp - Remove hairpin reads

## Synopsis
```
$ dnmtools cleanhp [OPTIONS] <read-1.fastq> <read-2.fastq>
```

## Description

## Options

```
 -o, -output
```
output filename prefix [required]
```
 -s, -stat
```
 stats output filename [required]
```
 -h, -hairpin
```
maximum hairpin rate
```
 -check
```
check for hairpin contamination
```
 -n, -nreads
```
 number of reads in initial check
```
 -c, -cutoff
```
 cutoff for calling an inverse duplication(default: 0.95)
```
 -i, -ignore
```
length of read name suffix to ignore when matching
```
 -v, -verbose
```
print more run info to STDERR while the program is running

