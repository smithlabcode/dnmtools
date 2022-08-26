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

The `uniq` command can take reads sorted by (chrom, start, end,
strand). If the reads in the input file are not sorted, run the
following sort command using [samtools](https://samtools.github.io):

```shell
$ samtools sort -O sam -o input-sorted.sam input.sam
```

Next, execute the following command to remove duplicate reads:

```shell
$ dnmtools uniq -S duplicate-removal-stats.txt input-sorted.sam out-sorted.sam
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
 -s, -seq
```
Use the sequences of the reads to distinguish duplicates. This is not
often recommended.

```txt
 -A, -all-cytosines
```
Use all cytosines when comparing reads based on sequence (default:
only use CpG sites). Only applies if `-s` (above) is used.

```txt
 -D, -disable
```
Disable testing if the reads are sorted by chromosome and
position. This can be faster and is fine if you know your reads are
sorted.

```txt
 -s, -seed
```
Random number seed. Which read to keep, among duplicates, is chosen
randomly (default: 408). This option is typically only used for
testing.

```txt
 -v, -verbose
```
Report more information while the program is running.
