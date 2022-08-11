# guessprotocol - Identify bisulfite sequencing protocol

## Synopsis
```
$ dnmtools guessprotocol [OPTIONS] <file-1.fastq> [<file-2.fastq>]
```

## Description

Mapping a WGBS dataset requires knowledge of the sequencing protocol
generated to process the data. This may not be properly documented
where the data was obtained, so we created a tool that guesses it
based on the nucleotide content on one or two FASTQ files.

The `guessprotocol` tool counts the number of As, Cs, Gs and Ts in
each end of the dataset and reports the protocol that is closest to
the nucleotide frequency expectations. outputs a single line
determining if the dataset is
 * (1) T-rich, where, in end 1, half of the bases are Ts and there are
very few Cs. In end 2, half of the bases are As and there are very few
Gs
 * (2) A-rich, where, in end 1, half of the bases are As and there
    are very few Gs and, in end 2, half of the bases are Ts and there
are very few Cs
 * (3) Random PBAT, where complementary reads have complementary
    bisulfite bases (e.g. if end 1 is T-rich, end 2 is A-rich), the
    bisulfite base in each end is random.
 * (4) Unknwon, if, based on the nucleotide frequencies, the protocol
   cannot be determined (or the read is not WGBS).

The output of `guessprotocol` is useful prior to mapping. For example,
it can be used to decide whether or not to map with the `-R` flag (for
"random PBAT") when using
[abismal](https://github.com/smithlabcode/abismal).

For paired-end data, `guessprotocol` finds read mates by finding
identical read names. Some datasets finish the read name with
identifiers like .1 on end 1 and .2 on end 2, thus making the read
names technically different at the last two characters. You can tell
the program to ignore a certain suffix size (like size 2 in this
example) when matching read names using the `-i` flag.

## Options
```
 -n -nreads
```
number of reads in initial check. The program stops after collecting
statistics for the first `n` reads (default: 1,000,000)

```
 -i -ignore
```
length of the read name suffix to ignore when matching
## Options
