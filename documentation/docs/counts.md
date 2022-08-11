# counts - compute single-site methylation

## Synopsis
```
$ dnmtools counts [OPTIONS] -c <chroms> <input.sam>
```


## Description

The `counts` command takes the mapped reads and produces the
methylation level at each genomic cytosine (both strands), with the
option to produce only levels for CpG-context cytosines.  While most
DNA methylation in vertebrates exists in the CpG context, cytosines in
other sequence contexts, such as CXG or CHH (where H denotes adenines,
thymines, or cytosines and X denotes adenines or thymines) may also be
methylated. Non-CpG methylation occurs most frequently in plant
genomes and pluripotent mammalian cells such as embryonic stem cells.
This type of methylation is asymmetric since the cytosines on the
complementary strand do not necessarily have the same methylation
status.

The input file (`input.sam`) is in SAM/BAM format. The reads should be
sorted so those mapping to the same chromosome are consecutive in the
file and [duplicate reads should be removed](../uniq).

The methylation level for every cytosine site at single base
resolution is estimated as the ratio of methylated to total bases in
reads mapped over that site. Because not all DNA methylation contexts
are symmetric, methylation levels are produced for both strands and
can be analyzed separately. To compute methylation levels at each
cytosine site along the genome you can use the following command:

```shell
$ dnmtools counts -c /path/to/genome.fa -o output.meth input.sam
```

The argument `-c` gives the filename of the genome sequence or the
directory that contains one FASTA format file for each chromosome. By
default `counts` identifies these chromosome files by the extension
`.fa`. Importantly, the "name" line in each chromosome FASTA file must
be the character `>` followed by the same name that identifies that
chromosome in the SAM output (the `.sam` files). An example of the
output and explanation of each column follows:

```txt
chr1    3000826 +   CpG 0.852941    34
chr1    3001006 +   CHH 0.681818    44
chr1    3001017 -   CpG 0.609756    41
chr1    3001276 +   CpGx    0.454545    22
chr1    3001628 -   CHH 0.419753    81
chr1    3003225 +   CpG 0.357143    14
chr1    3003338 +   CpG 0.673913    46
chr1    3003378 +   CpG 0.555556    27
chr1    3003581 -   CHG 0.641026    39
chr1    3003639 +   CpG 0.285714    7
```

The output file contains one line per cytosine site. The first column
is the chromosome. The second is the position of the cytosine
(zero-based). The 3rd column indicates the strand, which can be either
'+' or '-', corresponding to either a C on the positive reference
strand or a G on the positive reference strand. The 4th column is the
sequence context of that site, followed by an x if the site has
mutated in the sample away from the reference genome. The estimate for
mutation is based on whether the reads mapping on the opposite strand
for a cytosine are biased towards 'A', which happens if a 'C' to 'T'
mutation has occurred (a G to A on the other strand) and would
otherwise be interpreted as a C-to-T conversion by bisulfite. The 5th
column is the estimated methylation level, equal to the number of Cs
in reads at position corresponding to the site, divided by the sum of
the Cs and Ts mapping to that position. The final column is number of
reads overlapping with that site.

Note that because `counts` produces a file containing one line for
every cytosine in the genome, the file can get quite large. For
reference assembly mm10, the output is approximately 25GB. The `-n`
option produces methylation data for CpG context sites only, and for
mm10 this produces an output file that is approximately 1GB. It is
recommended that users allocate at least 8GB of memory when running
`counts`.

To examine the methylation status of cytosines a particular sequence
context, one may use the `awk` command (or `grep`) to filter those
lines based on the fourth column. For example, in order to pull out
all cytosines within the CHG context, run the following:

```shell
$ awk '$4 == CHG' human_esc.meth > human_esc_chg.meth
```

Our convention is to name `counts` output with all cytosines like
`*.meth`, with CHG like `*.CHG.meth` and with CHH like `*.CHH.meth`.

### Creating UCSC Genome Browser tracks

To view the methylation level or read coverage at individual CpG sites
in a genome browser, one needs to create a `bigWig` format file from a
counts file.

To create methylation level tracks or read coverage tracks, one can
follow these steps:

* (1) Download the wigToBigWig program from UCSC genome browser’s
  [directory of binary utilities](http://hgdownload.cse.ucsc.edu/admin/exe/).
* (2) Use the `fetchChromSizes` script from the same directory to
  create the \*.chrom.sizes file for the UCSC database you are
  working with (e.g. hg38). Note that this is the file that is
  referred to as `hg19.chrom.sizes` in step 3.
* (3) To create a `bigWig` track for methylation level at single CpG sites,
  convert the meth file to bed format using:
```shell
 $ awk -v OFS="\t" '{print $1, $2, $2+1, $4":"$6, $5, $3}' human_esc.meth > human_esc.meth.bed
```
* (4) To create a `bigWig` track from the bed format methcounts
  output, modify and use the following command:
```
 $ cut -f 1-3,5 human_esc.meth.bed | wigToBigWig /dev/stdin hg19.chrom.sizes human_esc.meth.bw
```
* (5) To create a `bigWig` track for coverage at single CpG sites, modify
  and use the following command:
```
$ tr ':' '[Ctrl+v Tab]' < human_esc.meth.bed | cut -f 1-3,5 | wigToBigWig /dev/stdin hg19.chrom.sizes human_esc.reads.bw
```

Note that if the `wigToBigWig` or `fetchChromSizes` programs are not
executable when downloaded, do the following:

```shell
$ chmod +x wigToBigWig
$ chmod +x fetchChromSizes
```

## Options

```txt
 -o, -output
```

Name of output file (default: stdout)

```txt
 -c, -chrom
```

File or directory of chroms (FASTA format; .fa suffix) [required].  If
the input is a directory, it should contain several FASTA files, each
one of which contains a chromosome sequence.

```txt
 -s, -suffix
```

Suffix of FASTA files (assumes -c specifies directory).

```txt
 -n, -cpg-only
```

Print only CpG context cytosines. This significantly reduces the
output size in most genomes.

```txt
 -v, -verbose
```

Report more information while the program is running.
