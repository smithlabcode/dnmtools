# guessprotocol - Identify bisulfite sequencing protocol

## Synopsis
```shell
$ dnmtools guessprotocol [OPTIONS] <file-1.fastq> [<file-2.fastq>]
```

## Description

Mapping a WGBS dataset requires knowledge of the sequencing protocol
generated to process the data. This may not be properly documented
where the data was obtained, so we created this command to guess the
protocol based on the nucleotide content in the input FASTQ file (or
files, for paired-end).

The `guessprotocol` tool uses two models of nucleotide content
following bisulfite conversion and applies this model to each
read. One model is for WGBS, and the other is for PBAT. For each read,
both models are applied, and the result is a probability for whether
the read (or read pair) was generated using WGBS or PBAT. Once the
requested number of reads is processed, the aggregate results for all
reads are used to guess whether the protocol used to generate the data
was WGBS, PBAT or rPBAT. The criteria are roughly as follows: if most
of the reads look like they are from WGBS, then we conclude WGBS.  If
most of the reads look like they are from PBAT, then we conclude
PBAT. If the result is more towards the middle, then we conclude
rPBAT.

More details: the number of As, Cs, Gs and Ts differs depending on
WGBS (traditional WGBS or MethylC-seq), PBAT -- post bisulfite adaptor
tagging, or rPBAT (random PBAT).

* For WGBS, a single-end sequenced read should be T-rich, and if the
  data is paired-end, read1 is T-rich and read2 is A-rich.
* For PBAT, a single-end sequenced read should be A-rich, and if the
  data is paired-end, read1 is A-rich and read2 is T-rich.
* For rPBAT, we have a random mix of the above situations. However, in
  practice it seems almost never to be 50% each.

In most cases, when the data is WGBS or PBAT, it is very obvious which
is the protocol used.

As of dnmtools v1.4.1, `guessprotocol` will always make a conclusion,
but includes a confidence level.

The output of `guessprotocol` is useful prior to mapping. For example,
it can be used to decide whether or not to map with the `-R` flag (for
"random PBAT") when using
[abismal](https://github.com/smithlabcode/abismal).

For paired-end data, `guessprotocol` finds ensures reads are mates by
finding identical read names. Some datasets finish the read name with
identifiers like ".1" on end 1 and ".2" on end 2, thus making the read
names technically different at the last two characters. You can tell
the program to ignore a certain suffix size (like size 2 in this
example) when matching read names using the `-i` flag.

The output includes the following values in a YAML format:
* `protocol`: this is the guessed protocol (wgbs, pbat or rpbat) based
  on the content of the reads.
* `confidence`: indicates the level of confidence in the guess for the
  protocol (values: low or high).
* `layout`: indicates whether the supplied reads were paired or
  single-ended.
* `n_reads_wgbs`: the average number of reads (for single-ended reads)
  or read pairs (for paired reads) where read1 is determined by the
  model to be T-rich.
* `n_reads`: the number of evaluated reads or read pairs.
* `wgbs_fraction`: the probability that a read (for single-ended
  reads) or the read1 of a read pair (for paired reads) is T-rich.

## Options
```
-n, -nreads
```
Number of reads to check. The program stops after collecting
statistics for the first `n` reads (default: 1,000,000). Fewer than
the default are usually sufficient, but increase this value if you
suspect reads at the start of the file might be problematic.

```txt
 -i -ignore
```
Length of the read name suffix to ignore when matching read names to
ensure mates are correctly synchronized when the data is paired-end.

```
-b, -bisulfite
```
Assumed bisulfite conversion rate for the models (default: 0.98).

```
-H, -human
```
Use human genome nucleotide frequencies. A good assumption for samples
from a mammal.

```
-o, -output
```
The output file name.

```
-v, -verbose
```
Report available information during the run.
