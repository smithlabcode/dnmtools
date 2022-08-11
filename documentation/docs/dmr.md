# dmr - Differentially methylated regions

## Synopsis
```shell
$ dnmtools dmr [OPTIONS] <input.methdiff> <input-a.hmr> <input-b.hmr> <output-a-lt-b.bed> <output-b-lt-a.bed>
```

## Description

With [differential methylation scores](../diff) and [HMRs](../hmr) for
both methylomes available, DM regions (DMRs) can be calculated with
the dmr program. This program uses HMR fragments to compute DMRs.

```shell
$ dnmtools dmr input.methdiff input-a.hmr input-b.hmr output-a-lt-b.bed output-b-lt-a.bed
```

The DMRs are output to files `output-a-lt-b.bed` and
`output-b-lt-a.bed`. The former file contains regions with lower
methylation in `input-a.meth` while the latter has regions with lower
methylation in `input-b.meth`.  Let’s take a look an example output
for `output-a-lt-b.bed`:

```txt
chr1    3539447 3540231 X:12  0 +
chr1    4384880 4385117 X:6   1 +
chr1    4488269 4488541 X:3   2 +
chr1    4603985 4604344 X:10  2 +
chr1    4760070 4760445 X:8   1 +
```

The first three columns are the genomic coordinates of DMRs. The fourth
column contains the number of CpG sites that the DMR spans, and the
fifth column contains the number of significantly differentially
methylated CpGs in the DMR. So, the first DMR spans 12 CpGs, but
contains no significantly differentially methylated sites, while the
second DMR spans 6 CpGs and contains just one significantly
differentially methylated CpG site.

We recommend filtering DMRs so that each one contains at least some
CpGs that are significantly differentially methylated. This can be
easily done with the `awk` utility, available on virtually all Linux and
Mac OS systems. For example, the following command puts all DMRs
spanning at least 10 CpGs and having at least 5 significantly
differentially methylated CpGs into a file `output-filtered.bed`.
```shell
$ awk -F "[:\t]" '$5 >= 10 && $6 >= 5' output.bed >output-filtered.bed
```

### Comparing two small groups of methylomes

To compare two small groups of methylomes, one should combine the
methylomes in each group and then compute DMRs for the resulting pair
of methylomes as described above. The methylomes can be combined using
the program [merge](../../utils/merge).

## Options

```txt
 -v, -verbose
```
print more run info to STDERR while the program is running.

```txt
 -c -cutoff
```
Significance cutoff (default: 0.05)

