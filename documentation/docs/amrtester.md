# amrtester - resolve epi-alleles

## Synopsis
```
$ dnmtools amrtester [OPTIONS] <input.bed> <input.epiread>
```

In addition to [amrfinder](../../amr/amrfinder), which uses a sliding
window, the `amrtester` program tests for allele-specific methylation
in a given set of genomic intervals. The program can be run like this:

```
$ dnmtools amrtester -o output.amr -c /path/to/genome.fa intervals.bed input.epiread
```

This program works very similarly to `amrfinder`, but does not have
options related to the sliding window. This program outputs a score
for each input interval, and when the likelihood ratio test is used,
the score is the p-value, which can easily be filtered later.

## Options

```
 -o, -output
```
output file name (default: STDOUT)
```
 -c, -chrom
```
FASTA file or directory of chromosomes containing FASTA files [required]
```
 -i, -itr
```
 max iterations
```
 -v, -verbose
```
 print more run info
```
 -P, -progress
```
print more run info to STDERR while the program is running.
```
 -b, -bic
```
use Bayesian Information Criterion (BIC) to compare models
