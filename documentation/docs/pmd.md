# pmd - Partially methylated domains

## Synopsis
```
$ dnmtools pmd [OPTIONS] <input.meth>
```

Huge genomic blocks with abnormal hypomethylation have been
extensively observed in human cancer methylomes and more recently in
extraembryonic tissues like the placenta. These domains are
characterized by enrichment in intergenic regions or Lamina associated
domains (LAD), which are usually hypermethylated in normal tissues.
Partially methylated domains (PMDs) are not homogeneously
hypomethylated as in the case of HMRs, and contain focal
hypermethylation at specific sites. PMDs are large domains with sizes
ranging from 10kb to over 1Mb. Hidden Markov Models can also identify
these larger domains. The program `pmd`  is provided for their
identification, and can be run as follows:
```
$ dnmtools pmd -i 1000 -o output.bed input.meth
```

The underlying Hidden Markov Model for PMD identification is very
similar to that of HMR identification, with a few key differences.
Because PMDs are large, megabase-scale regions of methylation loss
with individual CpGs exhibiting high methylation variance, each step
in the HMM is the weighted average methylation level for a genomic
region rather than a single CpG site. For samples with high sequencing
coverage and depth, a bin size of 1000bp suffices most of the time.

The PMD program has a built-in bin-size selection method that chooses
the smallest bin size (in 500bp increments) such that bins have enough
observations for their average methylation to be confidently
estimated. The bin size can be fixed by specifying `-b`. We recommend
defaulting to 1000 iterations to ensure the Baum-Welch training
procedure converges.

The sequence of genomic bins is segmented into hypermethylation and
partial-methylation domains, where the latter are the candidate PMDs.
Further processing of candidate PMDs includes trimming the two ends of
a domain to the first and last CpG positions, and merging candidates
that are "close" to each other. Currently, we are using twice the bin
size as the merging distance. Development in later versions of the
`pmd` program will include randomization procedures for choosing
merging distance.

The output of the program is a BED file formatted as follows:

```
chr1  4083054   4756012   PMD0:3.455556:278:4.584455:260  673  +
chr1  4846663   5430747   PMD1:3.842801:276:5.246181:166  580  +
chr1  5463102   5912049   PMD2:3.839765:217:4.191114:282  448  +
chr1  11366900  11716004  PMD3:4.141159:240:3.132796:411  349  +
chr1  12781853  13024981  PMD4:3.499296:241:2.168534:3    227  +
chr1  14741633  15210084  PMD5:5.306778:173:5.009876:147  467  +
chr1  18150329  18718393  PMD6:5.247428:179:3.984616:45   569  +
chr1  18789404  19197536  PMD7:3.511930:168:5.494317:268  408  +
chr1  29655012  29877666  PMD8:4.887023:232:1.406238:162  221  +
chr1  30028301  30464612  PMD9:5.101033:16:2.987482:118   434  +
```

The first three columns give the chromosome, start, and end positions
of the identified PMD. The fourth column has an arbitrarily assigned
name for the PMD (counting up from zero) and boundary quality
information: specifically the likelihood and certainty of the boundary
at the PMD start (1&2) and boundary at the PMD end (3&4).
These values can be used relative to other PMDs to assess
the "sharpness" of the PMD boundaries, and the confidence we have in
the boundary estimate based on the coverage of the sample. NaN
boundary scores represent the case where the joint likelihood and
confidence score of a boundary cannot be calculated because one side
of the boundary has either zero coverage or no CpGs: this occurs when
PMDs cut off due to deserts, for instance. The fifth column is the
number of bins in the PMD, which is analogous to the number of CpGs
segmented in the HMR program. PMDs should be called using information
from both strands, so the last column is a placeholder.

In general, the presence of a single HMR would not cause the program
to report a PMD in that region. However, in cases where a number of
HMRs are close to each other, such as the promoter HMRs in a gene
cluster, the pmd program might report a PMD covering those HMRs. Users
should be cautious with using such PMD calls in their further studies.
In addition, not all methylomes have PMDs, some initial visualization
or summary statistics can be of help in deciding whether to use `pmd`
program on the methylome of interest.

In addition, calling HMRs in samples with PMDs can be difficult: PMDs
can obscure the sites we are trying to identify by providing an
alternative foreground methylation state to the focused, very low
methylation typically at promoter regions. A good workaround for this
is to call PMDs first, and then call HMRs separately inside and
outside of PMDs (e.g. using [selectsites](../selectsites) and
using the output of PMDs as BED input). This ensures that the
foreground methylation state learned by the HMM in both types of
background is the focused hypomethylation at CpG islands and promoter
regions.

### PMD detection in multiple samples

The `pmd` program allows for multiple samples to be provided at once,
and it models each emission distribution separately. This behavior is
not yet well understood, but does yield an estimate of the "composite"
PMD state across the samples input. Because each emission distribution
is learned separately with no concept of a normally distributed error
on the parameters, this mode is less applicable with technical
replicates of the same tissue but potentially useful when looking
across many tissue types where substantial sample-to-sample biological
variability is expected.

### PMD detection in microarray data

Our recent addition of a bin size selection method, coupled with tuning
the desert size, can lead to a rough estimation of PMDs in human
Illumina MethylationEPIC microarray data. The `-a` option can be used to
specify that the input is all microarray data, and will result in the
emission distributions being modeled as Beta rather than
Beta-binomial. One should keep in mind that to use this option will
require some tuning of the desert size and bin sizes (we have seen
large desert sizes > 100kb and bin sizes around 20kb perform well) and
is best paired with matching WGBS data.

## Options

```
 -o, -out
```
output file (default: stdout)
```
 -d, -desert
```
maximum distance between bins with data in PMD
```
 -f, -fixedbin
```
Value of the fixed bin size
```
 -b, -bin
```
Starting bin size
```
 -a, -arraymode
```
Input is microarray data (e.g. from MethylationEPIC)
```
 -i, -itr
```
 max number of iterations

```
 -v, -verbose
```
print more run info to STDERR while the program is running
```
 -D, -debug
```
print debug info to STDERR while the program is runningo
```
 -P, -params-in
```
HMM parameter files for individual methylomes (separated with comma)
```
 -r, -posteriors-out
```
write out posterior probabilities in methcounts format
```
 -p, -params-out
```
write HMM parameters to this file
```
 -s, -seed
```
specify random seed value

