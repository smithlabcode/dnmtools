# uniq - ensure reads are not duplicates

## Synopsis
```shell
$ dnmtools uniq [OPTIONS] <input-sorted.sam> [out-sorted.sam]
```

## Description

The `uniq` command removes PCR duplicates. Before calculating
methylation level, you should now remove duplicate reads, which in
wgbs data are typically identified by their mapping to identical
genomic locations. These reads are most likely PCR clones rather than
representations of distinct DNA molecules. The command `uniq` remove
such duplicates. It collects duplicate reads and/or fragments that
have identical sequences and are mapped to the same genomic location
(same chromosome, same start and end positions, and same strand), and
chooses a random one to be the representative of the original DNA
sequence.

*Note* As of dnmtools v1.2.5, the option to use the sequence of reads
when deciding if two reads are duplicates has been removed. In the
context of analyzing bisulfite sequencing reads, this has the danger
of introducing bias in downstream analyses. Also, in the same version
the test for sorted order of reads cannot be disabled. Empirical tests
showed very little improvement to speed when disabling this test.

The `uniq` command can take reads sorted by (chrom, start, end,
strand). If the reads in the input file are not sorted, run the
following sort command using [samtools](https://samtools.github.io):

```shell
$ samtools sort -o reads_sorted.bam reads.bam
```

Next, execute the following command to remove duplicate reads:

```shell
$ dnmtools uniq -S duplicate-removal-stats.txt reads_sorted.bam reads_uniq.bam
```

## Options

```txt
 -S, -stats
```
Save statistics on duplication rates to this file. The statistics are not
reported unless a file is specified here.

```txt
 -hist
```
Output a histogram of duplication frequencies into the specified file
for library complexity analysis.

```txt
 -B, -bam
```
The output is in BAM format. This is an option to help prevent
accidentally writing BAM format to the terminal or through a pipe that
expects plain text, e.g., SAM.

```txt
 -stdout
```
Write the output to standard out. This is not done by default even
without an output file given, because of the danger of writing BAM to
the terminal or through a pipe unexpectedly. It is possible to write
BAM redirected or through a pipe, but the `-stdout` argument is
required.

```txt
 -s, -seed
```
Random number seed. Affects which read is kept among duplicates. The
default seed is 408. This option is typically only used for testing.

```txt
 -v, -verbose
```
Report more information while the program is running.
