# states - Allele-specific methylation file format

## Synopsis
```shell
$ dnmtools states [OPTIONS] <input.sam>
```

## Description

All programs that calculate statistics related to ASM must take the
linked states of CpG sites within reads into account. Using full read
sequences for this purpose is inefficient, so we defined an
intermediate format, "epiread," for this purpose. The `states` command
will convert a BAM or SAM file of mapped reads into a "states" file in
the format used by `amrfinder` and `amrtester`.

The epiread format consists of three columns. The first column is the
chromosome name for the mapped read, the second is the "index" of the
first CpG in the read. This index is a number (starting with 0) to
indicate which CpG sites in the chromosome corresponds to the first
CpG site in the read. These are not nucleotide positions in the
genome. The final column in the epiread format is the sequence of
methylation states within the read. This sequence of states is
composed of 3 possible letters: C if the corresponding letter at that
CpG site in the mapped read is a C, and similar for T. Within this
state sequence, letters in mapped reads at positions corresponding to
CpG sites that are neither C nor T are encoded as N. Aside from the
"N" this is effectively a binary encoding of methylation states.

Here is an example showing how some lines of an epiread format file might
look:
```txt
chr1    1460    CCCCCCCC
chr1    1460    CCC
chr1    1461    TCTTNNNNTTCT
chr1    1468    CCCC
chr1    1469    CCC
chr1    1469    CCCT
chr1    1469    CCC
chr1    1469    CCCCCCT
chr1    1469    CCC
chr1    1470    CCCC
chr1    1471    CCCNNNNNNTCCC
chr1    1472    CCC
```
Those epireads with the "N" in the middle correspond to paired-end
reads with ends that are joined. It is important to use these as one
fragment because linking methylation states within a fragment, over as
large a distance as possible, helps the inference methods within both
`amrfinder` and `amrtester`.

The following is an example of how to run the `states` command:
```shell
$ dnmtools states -c /path/to/genome.fa -o output.epiread input.sam
```

## Options

```txt
 -o, -output
```
The name of the output file.

```txt
 -c, -chrom
```
FASTA file of chromosomes containing FASTA files [required].

```txt
 -v, -verbose
```
Print more run info to STDERR while the program is running.
