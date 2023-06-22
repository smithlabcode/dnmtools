# radmerge - Merge CpGs to differentially methylated regions

## Synopsis
```shell
$ dnmtools radmerge [OPTIONS] <radmeth-input.bed>
```

## Description

After running [radmeth](../radmeth) followed by
[radadjust](../radadjust), it is possible to further join individually
differentially methylated CpGs into differentially methylated
regions. This can be achieved with the command

```shell
$ dnmtools radmerge -p 0.01 radmeth-input.bed > output-dmrs.bed
```

The current algorithm is conservative: it joins neighboring
differentially methylated sites with p-value below 0.01 (set by the -p
parameter). The output format is

```txt
chrom    start    end   dmr    num-sites   meth-diff
```

Above, `num-sites` and `meth-diff` are the number of significantly
differentially methylated CpGs in the DMR and the estimated
methylation difference, respectively. Example output might look like
this:

```txt
chr1     57315   57721  dmr     10      -0.498148
chr1     58263   59009  dmr     27      -0.521182
chr1    138522  139012  dmr     13      -0.443182
chr1    149284  149444  dmr      7      -0.430453
chr1    274339  275254  dmr     18      -0.520114
```

Note that in addition to being conservative, the work done by
`radmerge` is very simple, and does not consider genomic distance
between neighboring sites. It will merge consecutive significant sites
into one interval no matter how distant are those sites on a
chromosome.

## Options

```txt
 -o, -output
```
Output file (default: stdout).

```txt
 -p, -cutoff
```
P-value cutoff (default: 0.01).
