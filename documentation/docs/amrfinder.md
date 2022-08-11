# amrfinder - Compute allelically methylated regions (AMRs)

## Synopsis
```shell
$ dnmtools amrfinder [OPTIONS] <input.epiread>
```

## Description

The program `amrfinder` scans the genome using a sliding window to
identify AMRs. For a genomic interval, two statistical models are
fitted to the reads mapped, respectively. One model (single-allele
model) assumes the two alleles have the same methylation state, and
the other (two-allele model) represents different methylation states
for the two alleles.  Comparing the likelihood of the two models, the
interrogated genomic interval may be classified as an AMR.  The
following command shows an example to run the program amrfinder and
takes as input an epireads file generated from
[stats](../states).

```shell
$ dnmtools amrfinder -c /path/to/genome.fa -o output.amr input.epiread
```

There are several options for running amrfinder.

 * The `-b` switches from
using a likelihood ratio test to BIC as the criterion for calling an
AMR.

 * The `-i` option changes the number of iterations used in the EM
procedure when fitting the models.

 * The `-w` option changes the size of
the sliding window, which is in terms of CpGs. The default of 10 CpGs
per window has worked well for us.

 * The `-m` indicates the minimum
coverage per CpG site required for a window to be tested as an AMR.
The default requires 4 reads on average, and any lower will probably
lead to unreliable results.

 * The `-g` parameter is used to indicate the maximum distance between
   any two identified AMRS; AMRs are often fragmented, as coverage
fluctuates, and spacing between CpGs means their linkage cannot be
captured by the model.  if two are any closer than this value, they
are merged. The default is 1000, and it seems to work well in
practice, not joining things that appear as though they should be
distinct. In the current version of the program, at the end of the
procedure, any AMRs whose size in terms of base-pairs is less than
half the "gap" size are eliminated. This is a hack that has produced
excellent results, but will eventually be eliminated (hopefully soon).

 * The `-C` parameter specifies the critical value for keeping windows
   as AMRs, and is only useful when the likelihood ratio test is the
used; for BIC windows are retained if the BIC for the two-allele model
is less than that for the single-allele model.  amrfinder calculates a
false discovery rate to correct for multiple testing, and therefore
most p-values that pass the test will be significantly below the
critical value.

 * The `-h` option produces FDR-adjusted p-values according to a
   step-up procedure and then compares them directly to the given
critical value, which allows further use of the p-values without
multiple testing correction.

 * The `-f` omits multiple testing correction entirely by not applying
   a correction to the p-values or using a false discovery rate cutoff
to select AMRs.

## Options

```txt
 -o, -output
```
output file name (default: STDOUT)
```txt
 -c, -chrom
```
FASTA file or directory of chromosomes containing FASTA files [required]
```txt
 -i, -itr
```
max iterations
```txt
 -w, -window
```
size of sliding window (default: 10)
```txt
 -m, -min-cov
```
min coverage per CpG to test windows (default: 4)
```txt
 -g, -gap
```
min allowed gap, in bp, between AMRs (default: 1000)
```txt
 -C, -crit
```
critical p-value cutoff (default: 0.01)
```txt
 -f, -nofdr
```
omits FDR multiple testing correction
```txt
 -h, -pvals
```
 adjusts p-values using Hochberg step-up
```txt
 -b, -bic
```
use Bayesian Information Criterion (BIC) to compare models
```txt
 -v, -verbose
```
print more run info to STDERR while the program is running.
```txt
 -P, -progress
```
print progress info
