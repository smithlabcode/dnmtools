# cleanhp - Remove hairpin reads

## Synopsis
```shell
$ dnmtools cleanhp [OPTIONS] <read-1.fastq> <read-2.fastq>
```

## Description

## Options

```txt
 -o, -output
```
output filename prefix [required]
```txt
 -s, -stat
```
stats output filename [required]
```txt
 -h, -hairpin
```
maximum hairpin rate
```txt
 -check
```
check for hairpin contamination
```txt
 -n, -nreads
```
number of reads in initial check
```txt
 -c, -cutoff
```
cutoff for calling an inverse duplication(default: 0.95)
```txt
 -i, -ignore
```
length of read name suffix to ignore when matching
```txt
 -v, -verbose
```
print more run info to the terminal while the program is running
```txt
 -h, -hist
```
write a histogram of hairpin matches to this file
