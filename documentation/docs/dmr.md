# dmr - differentially methylated regions

## Synopsis
```console
$ dnmtools dmr [OPTIONS] <diffs> <hmrs-a> <hmrs-b> <dmr-a-lt-b> <dmr-b-lt-a>
```

## Description

When [differential methylation scores](../diff) and [HMRs](../hmr) for
two samples/methylomes available, differentially methylated regions
(DMRs) can be calculated using the `dmr` command. This command uses
the HMR fragments as candidate intervals, and indicates the number of
sites within each HMR fragment show significant differential
methylation in the appropriate direction (low vs. high, in both
directions). The "fragments" here are obtained by using the symmetric
difference of the given HMRs. The output gives this symmetric
difference in two separate files, and these may all be called
DMRs. However, the user should decide whether the fraction of sites
within each interval that show *significant* differential methylation
is sufficient to call any of these an DMR for a particular
analysis. So we expect users will filter these intervals (more on this
below).

Note 1: the `dmr` command will work -- meaning it should not crash --
if the two sets of provided intervals (the `input-a.hmr` and
`input-b.hmr` in the synopsis above) are not actually HMRs. As long as
they are in [BED](https://en.wikipedia.org/wiki/BED_(file_format))
format and non-overlapping within each file, the `dmr` command should
run. However, the intended use was with HMRs, and the interpretation
may be difficult of some other intervals are used.

Note 2: we refer to "CpG sites" here, but in some settings, for
example when studying Arabidopsis, non-CpG sites might be of
interest. The `dmr` command can also work in these settings.

The `dmr` command requires 5 files to be specified in a fixed
order. These include (1) the file of methylation differences as
produced by the `diff` command, (2) the HMRs for the first methylome,
(3) the HMRs for the second methylome, (4) the output file for DMRs
lower in the first methylome, and (5) the output file for DMRs lower
in the second methylome. The output file format is 6-column BED.

In the following example, the letters `lt` in the two output files
(the two BED files as the final arguments) indicate the direction of
differential methylation:

```console
$ dnmtools dmr input.diff input-a.hmr input-b.hmr dmr-a-lt-b.bed dmr-b-lt-a.bed
```

The DMRs are output to files `dmr-a-lt-b.bed` and
`dmr-b-lt-a.bed`. The former file contains regions with lower
methylation in the original counts output file for sample/methylome
"a", which might have been named `input-a.meth`. The latter has
regions with lower methylation in the counts output file for methylome
"b". One of these files, say `dmr-a-lt-b.bed` might look as follows:

```txt
chr1    3539447 3540231 X:12  0 +
chr1    4384880 4385117 X:6   1 +
chr1    4488269 4488541 X:3   2 +
chr1    4603985 4604344 X:10  2 +
chr1    4760070 4760445 X:8   1 +
```

The first three columns give genomic coordinates of the
"fragments". In this case, these would be intervals covered by
`input-a.hmr` and not bey `input-b.hmr`. The 4th column contains the
number of CpG sites that this DMR spans, preceded by the symbols "X:"
since this 4th column is a "name" in the bed format; we avoid simply
giving a number for this column. The 5th column contains the number of
significantly differentially methylated CpGs in this DMR where the
direction of the difference is lower in methylome "a" than in "b". So,
the first DMR spans 12 CpG sites, but contains no significantly
differentially methylated sites; the second DMR spans 6 CpGs and
contains just one significantly differentially methylated CpG site.

We recommend filtering DMRs so that each one contains a sufficient
number of CpG sites and meets some threshold for the fraction of CpG
sites that have significant differential methylation. This can be
easily done with `awk`, available on virtually all Linux and macOS
systems. For example, the following command filters to keep DMRs
spanning at least 10 CpGs and having at least 5 significantly
differentially methylated CpGs, storing them in a file named
`dmr-a-lt-b-filtered.bed`.

```console
$ awk -F "[:\t]" '$5 >= 10 && $6 >= 5' dmr-a-lt-b.bed > dmr-a-lt-b-filtered.bed
```

Above, the `-F` argument indicates possible field separator
characters, either a tab or the colon. If, for some reason, the tabs
in the file `dmr-a-lt-b.bed` have been converted to spaces, this would
break. If the fraction of significant CpGs is deemed more important
than their absolute number, for example at least 50% showing
significant differential methylation, the following command can be
used:

```console
$ awk -F "[:\t]" '$5 >= 10 && $6/$5 >= 0.5' dmr-a-lt-b.bed > dmr-a-lt-b-filtered.bed
```

**Warning:** the first input file, which is output from the
[diffs](../diffs) command, is directional. The direction determines
the order that the two HMR files must be specified. In the synopsis
above, the order of "hmrs-a" and "hmrs-b" must match the order that
methylomes/samples "a" and "b" were specified when running
[diffs](../diffs) to obtain the first of the input files. In a typical
application, if these are swapped we expect virtually no significant
sites in the DMRs (and the 5th column in the outputs would always be 0
or very close). So `diffs` would have been run as follows:

```console
$ dnmtools diff -o a-before-b.diff input-a.meth input-b.meth
```

### Comparing two small groups of methylomes

To compare two small groups of methylomes, one should combine the
methylomes (that is, the output from [counts](../counts)) within each
group and then compute DMRs for the resulting pair of methylomes using
the approach described above. The counts output files can be combined
using the program [merge](../merge). However, to take advantage of
replicates or experimental design, use the [radmeth](../radmeth)
command instead.

## Parameters

```txt
 -v, -verbose
```
Print more information while the command is running.

```txt
 -c -cutoff
```
Cutoff on p-values to define significant differences for individual
sites (default: 0.05).
