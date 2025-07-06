# counts-nano - compute single-site methylation from nanopore data

## Synopsis
```console
$ dnmtools counts-nano [OPTIONS] -c <chroms> <input.bam>
```

## Description

The `counts-nano` command introduced in v1.5.0 is designed specifically to
generate DNMTools [counts](../counts) format files from nanopore data called
for the `5mCG_5hmCG` modification. Currently this is only supported for
methylation and hydroxymethylation called at CpG sites.

More documentation will come as this tool evolves, but for now:

- Most behavior is very similar to what you will find from [counts](../counts).
- Mutation information is not estimated by `nano-counts`.
- Currently this only works for CpG sites and when the only modified sites are
  marked as `C+m?` or `C+h?` in the `MM` field of each BAM/SAM read record.
- The first 6 columns of the output are the same as explained in the
  [counts](../counts) format, except the fraction for the 5th column is both
  5mC and 5hmC. The 7th column is for 5hmC alone and the 8th is for 5mC alone.
- The methylation levels will not result in integer values when multiplied by
  the number of reads because probabilities on modifications are used, so
  methylation levels for each site are expected values (the best estimates we
  can make), and do not use arbitrary cutoffs.
- Several other commands in DNMTools have been modified to use this form of
  expected methylation level, and behave as previously for bisulfite
  sequencing data, but have updated behavior when the data is from
  nanopore. The user does not need to specify the technology used.
- Some commands need to use a `-relaxed` flag to work with the additional
  columns in the output from `counts-nano` compared with `counts`. For
  commands without this option, simply do `cut -f1-6` on the output of
  `counts-nano` to remove those.

## Options

```txt
-o, -output
```
Output file name. The default is to write output to the terminal,
which is not useful unless commands are piped.

```txt
-c, -chrom
```
Reference genome file, which must be in FASTA format. This is
required.

```txt
-t, -threads
```

The number of threads to use. This is only really helpful if the input is BAM
(not helpful for SAM), and the output is to be zipped (see `-z` below). These
threads are used to decompress BAM input and compress gzip output. If only one
of these conditions holds, using more threads can still help. Because most
computation in `counts-nano` is processing reads sequentially, using too many
threads will have decreasing returns.

```txt
-z, -zip
```

The output should be zipped (in gzip format). This is not deduced by the
filename, but specifying this argument should be accompanied by using a `.gz`
filename suffix for the output.

```txt
-n, -cpg-only
```

Print only CpG context cytosines. This significantly reduces the output size
in most genomes. Note that using this option does not merge data as symmetric
CpGs.

```txt
-sym
```

This will turn on `-n, -cpg-only` automatically and will output symmetric CpG
sites, with each level including all counts and methylation levels as a
(weighted) average of both strands.

```txt
-H, -header
```

Add a header to the output file to identify the reference genome. This will be
in the form of "comment" lines beginning with `#`. This is not required for most
downstream processing, but is used by commands that check for consistency with
a reference genome.

```txt
-v, -verbose
```

Report more information while the program is running.

```txt
-progress
```
Show progress while the program is running.
