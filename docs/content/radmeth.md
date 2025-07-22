## radmeth - Differential CpG methylation

## Synopsis

```console
dnmtools radmeth -factor case -o output.txt design_matrix.txt proportion_table.txt
```

## Description

The differential-methylation detection method described in this section is
based on the beta-binomial regression. We recommend using this method when
more than three replicates are available in each group. For rapid differential
methylation analysis the regression-based method should be run on a large
machine with many cores. Note that the actual processing time depends on the
coverage of each methylome, the number of sites analyzed, and the number of
methylomes in the dataset.

### Generating proportion tables

The first step in the differential methylation analysis is to assemble a
proportion table containing read proportions for all target
methylomes. Dnmtools includes a program [merge](../merge) to generate a
proportion table from the given collection of methylomes. Suppose that we want
to compare methylomes named `control-a.meth`, `control-b.meth`,
`control-c.meth` to the methylomes `case-a.meth`, `case-b.meth`,
`case-c.meth`.  The proportion table can be created with the following
command:

```console
$ dnmtools merge -t -radmeth \
     control-a.meth control-b.meth control-c.meth \
     case-a.meth case-b.meth case-c.meth > proportion-table.txt
```

This is what `proportion-table.txt` looks like:

```txt
control_a       control_b       control_c       case_a  case_b  case_c
chr1:108:+:CpG 9       6       10      8       1       1       2       2       2       1       14      1
chr1:114:+:CpG 17      7       10      0       14      3       5       1       9       1       7       1
chr1:160:+:CpG 12      8       10      5       17      4       15      14      13      6       4       4
chr1:309:+:CpG 1       1       1       0       17      12      12      8       2       1       19      8
```

As indicated in the header, this proportion table contains information about 6
methylomes: 3 controls and 3 cases.  Each row of the table contains
information about a CpG site and a proportion of reads mapping over this site
in each methylome. For example, the first row describes a cytosine within a
CpG site located on chromosome 1 at position 108. This site is present in 9
reads in the methylome control a and is methylated in 6 of them. Note that the
dnmtools `merge` command adds methylomes into the proportion table in the
order in which they are listed on the command line.

Note: the header conforms to a traditional "data frame" from R, which does not
include a column name for the row names, which are the first column.

### Design matrix

The next step is to specify the design matrix, which describes the structure
of the experiment. For our running example, the design matrix looks like this:

```txt
base   case
control_a      1       0
control_b      1       0
control_c      1       0
case_a 1       1
case_b 1       1
case_c 1       1
```

The design matrix shows that samples in this dataset are associated with two
factors: base and case. The first column corresponds to the base factor and
will always be present in the design matrix. Think of it as stating that all
samples have the same baseline mean methylation level. To distinguish cases
from controls we add another factor case (second column). The `1`s in this
column correspond to the samples which belong to the cases group. You can use
this design matrix as a template to create design matrices for two-group
comparisons involving arbitrary many methylomes.

After creating the proportion table and the design matrix, we are now ready to
start the differential methylation analysis.  It consists of (1) regression,
(2) combining significance, and (3) multiple testing adjustment steps.

### Regression

Suppose that the `proportion-table.txt` and `design-matrix.txt` are as
described above. The regression step is run with the command

```console
$ dnmtools radmeth -factor case -o output.txt design-matrix.txt proportion-table.txt
```

The `-factor` parameter specifies the factor with respect to which we want to
test for differential methylation. The test factor is 'case', meaning that we
are testing for differential methylation between cases and controls. In the
output file `output.txt`, the last four columns correspond to the total read
counts and methylated read counts of the case group and control group,
respectively. The p-value (5-th column) is given by the log-likelihood ratio
test.

Depending on the distribution of methylated and unmethylated CpGs, some
p-values may not be calculated. If there are no reads that cover a given CpG,
the value `NA_LOW_COV` will be printed instead of a p-value. If all CpGs are
fully methylated and fully methylated in both case and control, we say the CpG
has an "extreme distribution". In these cases, we cannot reject the null and
output `NA_EXTREME_CNT`.  Finally, if there are numerical errors in the
log-likelihood ratio test, we simply print `NA`. You can opt out of these
additional NA informations by not using the `-n` flag.

We do not recommend relying too much on the individual p-values of CpGs from
radmeth, as they often do not contain sufficient coverage for significant
biological interpretations. Instead, we recommend merging p-values of
neighboring CpGs using [radadjust](../radadjust).

## Tutorial

This tutorial aims to answer questions users often have about the assumptions
of radmeth and how to interepret the results. We will assume a factor of
interest, with levels A and B. We will also assume a potentially confounding
factor, sex, with levels M and F. Further, we have 8 total samples, 2 samples
for each of the 4 combinations of {A,B}x{M,F}. The design matrix should look
like this:

```txt
base    sex factor
sample_FA1  1   1   1
sample_FA2  1   1   1
sample_FB1  1   1   0
sample_FB2  1   1   0
sample_MA1  1   0   1
sample_MA2  1   0   1
sample_MB1  1   0   0
sample_MB2  1   0   0
```

Note that the 3 headings correspond with the three binary columns, but need
not be aligned over them. The first sample in the design matrix, `sample_FA1`,
has an indicator in the other columns for both the "F" and the "A", both
encoded with a "1" in their respective columns. But the replicate number 1 in
the `_FA1` is just a replicate, and if it had meaning we would consider
encoding it using some factor.  Designating it as a replicate with no meaning
is a choice here, aligning with the idea that the replicates are statistically
identical. If all replicates identified with a "2" were done on "day 2", for
example, then we might reconsider.

For our example we will use a genome file `genome.fa` that is very small and
not a real genome. We will also use a file `features.bed` that will be the
starting point for generating data that fit with our hypotheses. We will be
generating data that *should* show us the differences we want to see --
assuming we use the tools properly. I will refer to each of the intervals in
`features.bed` as a "feature" but just think of it as a genomic interval with
a chromosome name, a start and an end.  We need to first obtain some features
that are associated specifically with each of the levels of our factors: A, B,
M and F. We can do this starting with our features file and using the
`bedtools` shuffle command:

```console
$ bedtools shuffle -i features.bed -g genome.chrom.sizes > sim/features_M.bed;
$ cat sim/features_*.bed | sort -k 1,1 -k 2,2g -o excl.bed;
$ bedtools shuffle -excl excl.bed -i features.bed -g genome.chrom.sizes > sim/features_F.bed;
$ cat sim/features_*.bed | sort -k 1,1 -k 2,2g -o excl.bed;
$ bedtools shuffle -excl excl.bed -i features.bed -g genome.chrom.sizes > sim/features_A.bed;
$ cat sim/features_*.bed | sort -k 1,1 -k 2,2g -o excl.bed;
$ bedtools shuffle -excl excl.bed -i features.bed -g genome.chrom.sizes > sim/features_B.bed;
```

This ensures none of the "features" are overlapping between A, B, M or F, and
all are similar in number and size distribution. Using each of these, I then
made files of features corresponding to the *combinations* of factors:

```shell
for i in M F; do
    for j in A B; do
        cat sim/features_${i}.bed sim/features_${j}.bed |\
            sort -k 1,1 -k 2,2g -o sim/features_${i}${j}.bed;
    done;
done;
```

Now I have files named:

```txt
sim/features_MA.bed
sim/features_MB.bed
sim/features_FA.bed
sim/features_FB.bed
```

The files named with an "M" include all intervals in `sim/features_M.bed`, and
similarly the files named with a "B" include all intervals in
`sim/features_B.bed`. So the file `sim/features_FA.bed` is all the features
for "A" and all the features for "F". It will differ from
`sim/features_FB.bed` less than it will from `sim/features_MB.bed`.

Using these files, I simulated methylation levels at every CpG site in the
fake genome `genome.fa` with CpGs outside the intervals receiving a
methylation level of 0.8 and CpGs inside the intervals receiving a methylation
level of 0.15. At each CpG site, a number of reads was sampled (Poisson with
mean 10) and then the methylation levels were sampled with a coin flip weighed
as either 0.8 or 0.15, depending on the CpG site. I did this twice for each
combination of {M,F}x{A,B}.  The result is 8 files named as follows:

```console
$ ls -1 sim/*.counts
sim/sample_FA1.counts
sim/sample_FA2.counts
sim/sample_FB1.counts
sim/sample_FB2.counts
sim/sample_MA1.counts
sim/sample_MA2.counts
sim/sample_MB1.counts
sim/sample_MB2.counts
```

The program I used for the simulation above is not currently available. In
this simulation, the only difference between replicates (e.g.,
`sim/sample_FA1.counts` vs. `sim/sample_FA2.counts`) is random noise due to
sampling.

Here is a look inside one of these files:

```console
$ head sim/sample_FA1.counts
chr1    163 +   CG  0.857143    7
chr1    206 +   CG  0.928571    14
chr1    232 +   CG  0.7 10
chr1    278 +   CG  0.555556    9
chr1    296 +   CG  0.714286    7
chr1    310 +   CG  0.545455    11
chr1    322 +   CG  0.833333    6
chr1    324 +   CG  0.555556    9
chr1    350 +   CG  0.933333    15
chr1    356 +   CG  0.923077    13
```

Then I used the `merge` command from dnmtools to make a data matrix as
follows:

```console
$ dnmtools merge -t -radmeth -v -remove .counts -o sim/table.txt sim/sample_*.counts
```

And this is what it gave me:

```console
$ head sim/table.txt
sample_FA1  sample_FA2  sample_FB1  sample_FB2  sample_MA1  sample_MA2  sample_MB1  sample_MB2
chr1:163:+:CG   7   6   7   2   9   7   10  8   11  7   10  7   12  10  14  11
chr1:206:+:CG   14  13  4   2   12  11  15  14  8   7   11  8   5   4   11  9
chr1:232:+:CG   10  7   5   3   9   8   14  9   10  8   3   2   2   1   10  8
chr1:278:+:CG   9   5   16  14  14  12  8   8   10  9   8   7   6   4   9   7
chr1:296:+:CG   7   5   15  14  9   7   11  10  5   5   6   5   9   8   9   7
chr1:310:+:CG   11  6   11  11  8   6   9   5   5   4   10  7   9   7   5   4
chr1:322:+:CG   6   5   6   6   9   8   19  14  13  11  13  11  9   6   10  7
chr1:324:+:CG   9   5   10  9   8   5   10  9   7   6   14  9   7   6   6   6
chr1:350:+:CG   15  14  8   6   6   5   9   6   7   6   6   4   9   8   11  10
```

With this simulated data, in which each column has the influence of two
different factors, we can use `radmeth` with the design matrix shown above (in
a file named `desgin.txt`) to get p-values for differential methylation
between levels A and B of our factor of interest -- while controling for the
effect of M/F:

```console
$ dnmtools radmeth -o sim/samples.radmeth -f factor design.txt sim/table.txt
```

As explained already, the 5th column of the output gives the p-values for
differential methylation between levels for our factor of interest (named
"factor" in the design matrix).  Here is a look at part of the output:

```console
$ head sim/samples.radmeth
chr1    163 +   CG  0.0864953   35  22  45  36
chr1    206 +   CG  0.444024    37  30  43  38
chr1    232 +   CG  0.756591    28  20  35  26
chr1    278 +   CG  0.780385    43  35  37  31
chr1    296 +   CG  0.64704 33  29  38  32
chr1    310 +   CG  0.640826    37  28  31  22
chr1    322 +   CG  0.0971364   38  33  47  35
chr1    324 +   CG  0.257851    40  29  31  26
chr1    350 +   CG  0.893476    36  30  35  29
chr1    356 +   CG  0.0527032   37  35  42  33
```

Focusing on the 5th column, we see values between 0 and 1, but none of them
are very extreme. This particular output file has 17903 rows.  We can use the
dnmtools command `selectsites` to pull out the rows in this file that
correspond to each of the original feature sets. Here is how we would do it
for all 4 different original sets of features:

```shell
for i in M F A B; do
    dnmtools selectsites -o sim/features_${i}.txt sim/features_${i}.bed sim/samples.radmeth;
done
```

The result allows us to verify that `radmeth` with the design matrix given
above identifies sites as significant when they differ between levels A and B
of the factor of interest:

```console
$ for i in sim/features_*.counts; do echo $i `awk 'BEGIN{k=0;c=0}{k+=$5;c+=1}END{print k,c,k/c}' $i`; done | column -t
sim/features_A.counts  0.255448  430  0.000594064
sim/features_B.counts  0.114256  372  0.00030714
sim/features_F.counts  196.093   400  0.490232
sim/features_M.counts  295.234   588  0.502099
```

The p-values in features associated with either M or F are averaging around
0.5, which is what we want if we are trying to control this possible
confounding variable. In the features associated either with A or B, the
p-values are very small on average. Of course this is literally by
design. Knowing that we simulated just as much difference between M and F as
between A and B, if we had used `-f sex` then `radmeth` would have correctly
assigned low p-values to sites in the intervals of `features_M.bed` and
`features_F.bed`.

We can use the same simulated data to verify that we do need the intercept in
our design. In other words, the column "base" in `design.txt` needs to be
present (you can give it any name you want; make sure it's all 1s).

We change the design matrix to make the following, which has the intercept
column removed:

```console
$ cat design_noint.txt
sex factor
sample_FA1  1   1
sample_FA2  1   1
sample_FB1  1   0
sample_FB2  1   0
sample_MA1  0   1
sample_MA2  0   1
sample_MB1  0   0
sample_MB2  0   0
```

Then we re-run all the above commands (I won't exaplain each) to get results
analogous but with this modified design matrix:

```console
$ dnmtools radmeth -o sim/samples_noint.radmeth -f factor design_noint.txt sim/table.txt
$ for i in M F A B; do dnmtools selectsites -o sim/features_noint_${i}.txt sim/features_${i}.bed sim/samples_noint.radmeth; done
$ for i in sim/*_noint_*.txt; do echo $i `awk 'BEGIN{k=0;c=0}{k+=$5;c+=1}END{print k,c,k/c}' $i`; done | column -t
sim/features_noint_A.txt  8.27788  430  0.0192509
sim/features_noint_B.txt  15.2391  372  0.0409653
sim/features_noint_F.txt  77.2438  400  0.193109
sim/features_noint_M.txt  93.7181  588  0.159385
```

Although the sites within features corresponding to A and B are much lower on
average than those for M or F, the distinction has been reduced. More
importantly, we know that the sites in `sim/features_M.bed` and
`sim/features_F.bed` are not influenced by the level A/B of our factor of
interest. The p-values for sites in the last two rows above should have mean
of 0.5. But without the intercept included in the design matrix, the outcomes
to not match our statistical assumptions. I cannot predict what will happen if
you use the wrong design matrix. In analyzing bisulfite sequencing data for
differential methylation, we become aware of uncertainty through the counts,
and we leverage this to try and be more accurate. So unlike other regression
settings, if you try to normalize the data first, in order to avoid having to
deal with an intercept, you will not arrive at an equivalent procedure -- if
you want to normalize the data yourself (i.e. subtract mu and divide by sigma)
then you are probably better off using a general regression package.

## Options

```txt
 -o, -out
```

The name of the output file (default: stdout).

```txt
 -n -na-info
```

If a p-value is not calculated, print NAs in more detail: low coverage
(`NA_LOW_COV`), extreme values (`NA_EXTREME`) or numerical errors in
likelihood ratios (`NA`).

```txt
 -v, -verbose
```

Print more information while the command is running.

```txt
 -progress
```

Show progress while radmeth runs.

```txt
 -f, -factor
```

The name of the factor on which to test differences. This must be the
name of one of the columns in file design matrix. This is required.

```txt
 -t, -threads
```

Use multiple threads. This gives very good speedup as long as you do not
exceed the number of available physical cores. Most CPUs support 2 hardware
threads per physical core via SMT/Hyper-Threading, so check your number of
cores if you are unsure.
