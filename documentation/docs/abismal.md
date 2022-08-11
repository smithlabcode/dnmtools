# abismal - another bisulfite mapping algorithm

## Synopsis
```shell
$ dnmtools abismal [OPTIONS] input.fq [input-r2.fq]
```

# Description
During bisulfite treatment, unmethylated cytosines in the original DNA
sequences are converted to uracils, which are then incorporated as
thymines (T) during PCR amplification. These PCR products are referred
to as T-rich sequences as a result of their high thymine constitution.
With paired-end sequencing experiments, the compliments of these
T-rich sequences are also sequenced.  These complimentary sequences
have high adenine (A) constitution (A is the complimentary base pair
of T), and are referred to as A-rich sequences. Mapping consists of
finding sequence similarity, based on context specific criteria,
between these short sequences, or reads, and an orthologous reference
genome.  When mapping T-rich reads to the reference genome, either a
cytosine (C) or a thymine (T) in a read is considered a valid match
for a cytosine in the reference genome. For A-rich reads, an adenine
or a guanine is considered a valid match for a guanine in the
reference genome. The mapping of reads to the reference genome by
`abismal` is described below. If you choose to map reads with a
different tool, make sure that your post-mapping files are
appropriately formatted for the next components of the `dnmtools`
pipeline (necessary file formats for each step are covered in the
corresponding sections).  The default behavior of `abismal` is to
assume that reads are T-rich and map accordingly, but different
sequencing protocols that generate A-rich and T-rich reads in
different combinations are equally accepted. `abismal` is
available [here](http://github.com/smithlabcode/abismal).

### Input and output file formats

We assume that the original data is a set of sequenced read files,
typically as produced by Illumina sequencing. These are FASTQ format
files, and can be quite large. After the reads are mapped, these files
are not used by our pipeline. The reference genome should be a single
FASTA file that was previously indexed using the `abismalidx`
tool. The `abismal` program requires an indexed FASTA reference genome
and the input FASTQ files(s), after which it generates a Sequence
Alignment/Map(SAM) output indicating the coordinates of mapped reads.
Details of the SAM file format can be found at the SAM file format
[documentation](http://samtools.github.io/hts-specs/SAMv1.pdf). These
SAM files will be the input files for the postprocessing quality
control and analysis programs to follow, including [bsrate](../bsrate)
and [counts](../counts).

`abismal` operates by preprocessing the reference genome into a large
index, where k-mers of set length are used as keys to a list of
potential mapping locations for reads that begin with their sequence.
To produce this index run the following command:

```shell
$ abismalidx  <genome folder or file> <index file>
```

For the human genome, whose size is 3 GB, the resulting index is
approximately 2.5 GB.

### Decompressing and isolating paired-end reads

Sometimes paired-end reads are stored in the same FASTQ file.  Because
we treat these paired ends differently, they must be separated into
two files and both files must be passed as inputs to `abismal`.

If your data is compressed as a Sequenced Read Archive, or SRA file,
you can decompress and split paired-end reads into two files at the
same time using `fastq-dump` , which is a program included in the
[sra-toolkit](https://hpc.nih.gov/apps/sratoolkit.html)
package, available for most unix systems.  Below is an example of using
`fastq-dump` to decompress and separate FASTQ data by end:

```shell
$ fastq-dump --split-3 human_esc.sra
```

If you have a FASTQ file not compressed in SRA format, you can split
paired ends into two separate files by running the following commands:

```shell
$ sed -ne '1~8{N;N;N;p}' *.fastq > *_1.fastq
$ sed -ne '4~8{N;N;N;p}' *.fastq > *_2.fastq
```

### Sequencing adapters

These are a problem in any sequencing experiment with short fragments
relative to the lengths of reads. A robust method for removing
adapters is available via the Babraham Institute's Bioinformatics
group, called
[trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore).
This program takes into account adapter contamination of fewer than 10
nucleotides and can handle mismatches with the provided sequence. We
strongly recommend reads are trimmed prior to mapping. Our recommended
parameter choice for `trim_galore` restricts trimming only to
sequencing adapters, leaving all other bases unaltered. This can be
attained by running it with the following command for single-end
mapping

```shell
$ trim_galore -q 0 --length 0 human_esc.fastq
```

and the following command for paired-end

```shell
$ trim_galore --paired -q 0 --length 0 human_esc_1.fastq human_esc_2.fastq
```

Note that the `--paired` flag is necessary to ensure the program does
not interpret the two input files as independent and that the
resulting FASTQ files still have corresponding read mates in the
correct order.

### Single-end reads

When working with data from a single-end sequencing experiment, you
will have T-rich reads only. `abismal` expects T-rich reads as a
default. Execute the following command to map all of your single-end
reads with `abismal`:

```shell
$ abismal -i <genome index> -o <output SAM> [options] <input fastq>
```

To save files in BAM format, which significantly reduce disk space,
simply redirect the `abismal` output to the [samtools
view](https://github.com/samtools/samtools) program using the `-b`
flag to compress to BAM and the `-S` flag to indicate that input is in
SAM format.

```shell
$ abismal -i <genome index> -s <output STATS> [options] <input fastq> | samtools view -bS > <output BAM>
```

### Paired-end reads

When working with data from a paired-end sequencing experiment, you
will have T-rich and A-rich reads. T-rich reads are often kept in
files labeled with an `_1` and A-rich reads are often kept in files
labeled with an `_2`.  T-rich reads are sometimes referred to as
5' reads or mate 1 and A-rich reads are sometimes referred
to 3' reads or mate 2. We assume that the T-rich file and
the A-rich contain the same number of reads, and each pair of mates
occupy the same lines in their respective files. We will follow this
convention throughout the manual and strongly suggest that you do the
same. Run the following command to map two
reads files from a paired-end sequencing experiment:

```shell
$ abismal -i <index> -o <output SAM> [options] <input fastq 1> <input fastq 2>
```

In brief, what happens internally in `abismal` is as follows.
`abismal` finds candidate mapping locations for a T-rich mate with
CG-wildcard mapping, and candidate mapping locations for the
corresponding A-rich mate with AG-wildcard mapping. If two candidate
mapping locations of the pair of mates are within certain distance in
the same chromosome and strand and with correct orientation, the two
mates are combined into a single read (after reverse complement of the
A-rich mate), referred to as a fragment. The overlapping region
between the two mates, if any, is included once, and the gap region
between them, if any, is filled with Ns. The parameters `-l` and `-L`
to `abismal` indicate the minimum and maximum size of fragments to
allow to be merged, respectively.

Here the fragment size is the sum of the read lengths at
both ends, plus whatever distance is between them. So this is the
length of the original molecule that was sequenced, excluding the
sequencing adapters. It is possible for a given read pair that the
molecule was shorter than twice the length of the reads, in which case
the ends of the mates will overlap, and so in the merged fragment will
only be included once. Also, it is possible that the entire molecule
was shorter than the length of even one of the mates, in which case
the merged fragment will be shorter than either of the read ends. If
the two mates cannot be merged because they are mapped to different
chromosomes or different strand, or they are far away from each other,
`abismal` will output each mate individually if its mapping
position is unambiguous.

`abismal` provides a statistical summary
of its mapping job in the `mapstats` file, which includes the total
number and proportion of reads mapped, how many paired end mates were
mapped together, and the distribution of fragment lengths computed by
the matched pairs. The `mapstats` file name is usually the same as
the SAM output with `mapstats` appended to it, but custom file
names can be provided using the `-s` flag.

### Mapping reads in a large file:

Mapping reads may take a while, and mapping reads from WGBS takes
even longer. It usually take quite a long time to map reads from a
single large file with tens of millions of reads. If you have access
to a cluster, one strategy is to launch multiple jobs, each working on
a subset of reads simultaneously, and finally combine their output.
`abismal` takes advantage of OpenMP to parallelize the process of
mapping reads using the shared index.

If each node can only access its local storage, dividing the set of
reads to into k equal sized smaller reads files, and mapping these all
simultaneously on multiple nodes, will make the mapping finish about k
times faster.  The UNIX `split` command is good for dividing the reads
into smaller parts. The following BASH commands will take a directory
named `reads` containing Illumina sequenced reads files, and split
them into files containing at most 3M reads:

```shell
$ mkdir reads_split
$ for i in reads/*.txt; do split -a 3 -d -l 12000000 ${i} reads_split/$(basename $i); done
```

Notice that the number of lines per split file is 12M, since we want
3M reads, and there are 4 lines per read. If you split the reads like
this, you will need to ``unsplit'' them after the mapping is done. This
can be done using the `samtools merge` command.

Abismal also exists as a standalon program, and more details on
abismal are available in its [documentation
manual](https://github.com/smithlabcode/abismal/blob/master/docs/MANUAL.md#quick-installation)
