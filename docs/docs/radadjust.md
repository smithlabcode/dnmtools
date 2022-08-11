# radadjust - Correct p-values of individual CpGs

## Synopsis
```
$ dnmtools radadjust [OPTIONS] <regression-output.bed>
```

## Description

This program adjusts the p-value of individual CpGs in the output of
[radmeth](../../dmr/radmeth). A typical application
that takes the regression output as input and combines the p-values of
200 neighboring CpGs is done as follows.
```
$ dnmtools radadjust -bins 1:200:1 input.bed >output-adjusted.bed
```

Here, the only required parameter, besides the input file, is `-bins`
whose value is set to `1:200:1` (which is also the default value). This
means that for each `n = 1, 2, ...199`, `radmeth-adjust` computes the
correlation between p-values of CpGs located at distance n from each
other. These correlations are used during significance combination
step. In addition, bin sizes determine the window for combining
significance. In contrast, if `-bins` is set to `1:15:5`, then the
correlation is computed separately for p-values corresponding to CpGs
at distances `[1, 5)`, `[5, 10)`, and `[10, 15)` from one another.  The
first five columns and the last four columns of `radmeth-adjust` have
the same meaning as those output by radmeth regression. The 6th column
gives the modified p-value based on the original p-value of the site
and the p-values of its neighbors. The 7th column gives the
FDR-corrected p-value. Then the last four columns correspond to the
total read counts and methylated read counts of the case group and
control group, respectively.  Here is what the `output-adjusted.bed`
file looks like for our example dataset:

```
chr1  108   +     CpG   0.157971    0.099290    0.353466    18     4    20    15
chr1  114   +     CpG   0.559191    0.099290    0.353466    21     3    41    10
chr1  160   +     CpG   0.095112    0.099290    0.353466    32    24    39    17
chr1  309   +     CpG   0.239772    0.122248    0.368902    33    17    19    13
chr1  499   +     CpG   0.770140    0.204467    0.419872    43    22    29    15
```

After completing the previous steps, individual differentially methy-
lated sites can be obtained with 'awk'. To get all CpGs with
FDR-corrected p-value below 0.01, run

```
$ awk '$7 <= 0.01' output-adjusted.bed >output-significant.bed
```

## Options

```
 -o, -out
```
output file (default: STDOUT)

```
 -b, -bins
```
correlation bin specs

```
 -v, -verbose
```
print more run info to STDERR while the program is running.

