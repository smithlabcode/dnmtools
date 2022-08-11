## radmeth - Differential CpG methylation

## Synopsis
```
$ dnmtools radmeth -factor case design_matrix.txt proportion_table.txt >output.bed
```

## Description

The DM detection method described in this section is based on the
beta-binomial regression. We recommend using this method when more
than three replicates are available in each group. For rapid
differential methylation analysis the regression-based method should
be run on a computing cluster with a few hundred available nodes, in
which case it takes approximately a few hours to process a dataset
consisting of 30-50 WGBS samples. The analysis can be also performed
on a personal workstation, but it will take substantially longer. Note
that the actual processing time depends on the coverage of each
methylome, the number of sites analyzed, and the number of methylomes
in the dataset.

### Generating proportion tables

The first step in the differential methylation analysis is to assemble
a proportion table containing read proportions for all target
methylomes. Dnmtools includes a program
[merge](../merge) to generate a proportion
table from the given collection of methylomes. Suppose that we want to
compare methylomes named `control-a.meth`, `control-b.meth`, `control-
c.meth` to the methylomes `case-a.meth, `case-b.meth, `case-c.meth.
The proportion table can be created with the following command:

```
$ dnmtools merge -t control-a.meth control-b.meth control-c.meth \

case-a.meth case-b.meth case-c.meth > proportion-table.txt
```

This is what `proportion-table.txt` looks like:

```
	control_a	control_b	control_c	case_a	case_b	case_c
chr1:108:+:CpG	9	6	10	8	1	1	2	2	2	1	14	1
chr1:114:+:CpG	17	7	10	0	14	3	5	1	9	1	7	1
chr1:160:+:CpG	12	8	10	5	17	4	15	14	13	6	4	4
chr1:309:+:CpG	1	1	1	0	17	12	12	8	2	1	19	8
```

As indicated in the header, this proportion table contains information
about 6 methylomes: 3 controls and 3 cases.  Each row of the table
contains information about a CpG site and a proportion of reads
mapping over this site in each methylome. For example, the first row
describes a cytosine within a CpG site located on chromosome 1 at
position 108. This site is present in 9 reads in the methylome control
a and is methylated in 6 of them. Note that `merge-methcounts` adds
methylomes into the proportion table in the order in which they are
listed on the command line.

### Design matrix

The next step is to specify the design matrix, which describes the
structure of the experiment. For our running example, the design
matrix looks like this:

```
base	case
control_a	1	0
control_b	1	0
control_c	1	0
case_a	1	1
case_b	1	1
case_c	1	1
```

The design matrix shows that samples in this dataset are associated
with two factors: base and case. The first column corresponds to the
base factor and will always be present in the design matrix. Think of
it as stating that all samples have the same baseline mean methylation
level. To distinguish cases from controls we add another factor case
(second column). The `1`s in this column correspond to the samples
which belong to the cases group. You can use this design matrix as a
template to create design matrices for two-group comparisons involving
arbitrary many methylomes.

After creating the proportion table and the design matrix, we are now
ready to start the methylation analysis.  It consists of (1)
regression, (2) combining significance, and (3) multiple testing
adjustment steps.

### Regression

Suppose that the `proportion-table.txt` and `design-matrix.txt` are as
described above. The regression step is run with the command

```
$ dnmtools radmeth -factor case design-matrix.txt proportion-table.txt > output.bed
```

The `-factor` parameter specifies the factor with respect to which we
want to test for differential methylation. The test factor is case,
meaning that we are testing for differential methylation between cases
and controls. In the output file `output.bed`, the last four columns
correspond to the total read counts and methylated read counts of the
case group and control group, respectively. The p-value (5-th column)
is given by the log-likelihood ratio test.

Depending on the distribution of methylated and unmethylated CpGs,
some p-values may not be calculated. If there are no reads that cover
a given CpG, the value `NA_LOW_COV` will be printed instead of a
p-value. If all CpGs are fully methylated and fully methylated in both
case and control, we say the CpG has an "extreme distribution". In
these cases, we cannot reject the null and output `NA_EXTREME_CNT`.
Finally, if there are numerical errors in the log-likelihood ratio
test, we simply print `NA`. You can opt out of these additional NA
informations by not using the `-n` flag.

We do not recommend using the individual p-values of CpGs from
radmeth, as they often do not contain sufficient coverage for
significant biological interpretations. Instead, we recommend merging
p-values of neighboring CpGs using [radadjust](../radadjust).

## Options

```
 -o, -out
```
output file (default: STDOUT)

```
 -n -na-info
```
if a p-value is not calculated, print NAs in more
detail: low count (`NA_LOW_COV`) extreme values (`NA_EXTREME`)
or numerical errors in likelihood ratios (`NA`).

```
 -v, -verbose
```
print more run info to STDERR as the program runs.

```
 -f, -factor
```
a factor to test, one of the columns in file `design-matrix.txt` [required]


