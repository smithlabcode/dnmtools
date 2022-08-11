# hmr-rep - Hypomethylated regions across replicates

## Synopsis
```shell
$ dnmtools hmr-rep [OPTIONS] <input-1.meth> <input-2.meth> ...
```

## Description

This program is similar to [hmr](../hmr), but it identifies
HMRs in a set of replicate methylomes. Methylation must be provided in
the [counts](../counts) format. This program assumes
only data at CpG sites and that strands are collapsed so only the
positive site appears in the file (e.g. using
[sym](../sym)).

## Options

```txt
 -o, -out
```
output file (default: STDOUT)
```txt
 -d, -desert
```
maximum distance between covered CpGs in HMR (default: 1000)

```txt
 -i, -itr
```
max number of iterations (default: 100)
```txt
 -v, -verbose
```
print more run info to STDERR while the program is running.
```txt
 -partial
```
identify PMRs instead of HMRs
```txt
 -post-hypo
```
output file for single-CpG posterior hypomethylation probability (default: none)

```txt
 -post-meth
```
output file for single-CpG posteiror methylation probability (default: none)

```txt
 -P, -params-in
```
HMM parameter file (override training step)
```txt
 -p, -params-out
```
write HMM parameters to this file (default: none)
```txt
 -s, -seed
```
specify random seed (default: 408)

