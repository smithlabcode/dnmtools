# counts - compute single-site methylation

## Synopsis
```console
$ dnmtools counts [OPTIONS] -c <chroms> <input.bam>
```

## Description

The `counts` command takes the mapped reads and produces the
methylation level at each genomic cytosine (both strands), with the
option to produce only levels for CpG-context cytosines.  While most
DNA methylation in vertebrates is in the CpG context, cytosines in
other sequence contexts, such as CXG or CHH (where H denotes adenines,
thymines, or cytosines and X denotes adenines or thymines) may also be
methylated. Non-CpG methylation occurs frequently in plant genomes and
pluripotent mammalian cells such as embryonic stem cells. And possibly
whatever cells you are studying. The output of `counts` serves as the
input for many downstream analyses.

The input mapped reads file (`input.bam`) is in SAM/BAM format. The
reads should be sorted so those mapping to the same chromosome are
consecutive in the file. Duplicate reads should be probably be
[removed](../uniq) first, but that depends on your data.

The methylation level for every cytosine site at single base
resolution is estimated as the ratio of methylated to total bases in
reads mapped over that site. If a site has no reads covering it, then
we leave the value at 0, but also indicate this in the output (see
below).

Each site is marked with a *context* in the output of the `counts`
command, and you can find out more about cytosine contexts
[here](../cytosine_contexts).

To compute methylation levels at each cytosine site along
the genome you can use the following command:
```console
$ dnmtools counts -c /path/to/genome.fa -o output.meth input.sam
```

The argument `-c` gives the name of a FASTA file containing all
chromosome sequences (as of v1.2.5, a directory of separate files is
no longer supported). Importantly, the "name" line in each chromosome
FASTA file must begin with the character `>` followed immediately by
the same name that identifies that chromosome in the SAM output (the
`.bam` files). If you use the same FASTA format file you used to map
the reads, everything should be fine. An example of the output and
explanation of each column follows:
```txt
chr1  1869  +  CCG  0         1
chr1  1870  +  CpG  0         1
chr1  1871  -  CpG  0.666667  6
chr1  1874  +  CHH  0         1
chr1  1875  +  CCG  0         1
chr1  1876  +  CpG  0         1
chr1  1877  -  CpG  0.333333  6
chr1  1880  +  CHH  0         1
chr1  1881  +  CCG  0         1
chr1  1882  +  CpG  0         1
chr1  1883  -  CpG  0.5       6
chr1  1886  +  CHH  0         1
chr1  1887  +  CCG  0         1
chr1  1888  +  CpG  0         1
chr1  1889  -  CpG  0.333333  6
chr1  1892  +  CHH  0         1
chr1  1893  +  CCG  0         1
chr1  1894  +  CpG  0         1
chr1  1895  -  CpG  0.25      4
chr1  1898  +  CHH  0         1
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
context, you can use the `awk` command (or `grep`) to filter those
lines based on the fourth column. For example, to get all cytosines
within the CHG context, run the following:
```console
$ awk '$4 == CHG' human_esc.meth > human_esc_chg.meth
```
Our convention is to name `counts` output with all cytosines like
`*.meth`, with CHG like `*.chg.meth` and with CHH like `*.chh.meth`.

**Estimated "mutations" in the sample:**
In version 1.2.1 sites are flagged if there is an indication that a
mutation at a C in the sample might cause a T to appear in the reads,
and then confound the issues of mutations vs. the absence of
methylation. The 4th column in the output (as indicated above) is the
"context" string for that cytosine, and if a site is flagged as a
possible mutation, then an `x` will appear at the end of the context
string. This is an example of what you would see in the output:

```txt
chr1  1870  +  CpGx  0         1
```

The data used to make this inference is on the opposite strand: if
there is a C on the positive strand, then we should observe a G at the
same position in reads mapping to the negative strand, regardless of
the bisulfite. If the reference has a C, but molecules in the sample
have mutated C->T, independent of methylation or bisulfite, then we
would observe an A in reads mapping to the same position on the
opposite strand. At present we use very poor statistical criteria for
this, so unless you believe you understand it well, it's best to
ignore. Hopefully we can do this in a more rigorous and interpretable
way in the future.

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
The number of threads to use. This is only really helpful if the input
is BAM (not very helpful for SAM), and the output is to be zipped (see
`-z` below). These threads will help decompress the BAM input and will
help compress the gzip format output. If only one of these conditions
holds, using more threads can still help. Because `counts` spends most
of its computing time processing reads sequentially, there are
diminishing returns for specifying too many threads.

```txt
-z, -zip
```
The output should be zipped (in gzip format). This is not deduced by
the filename, but specifying this argument should be accompanied by
using a `.gz` filename suffix for the output.

```txt
-n, -cpg-only
```
Print only CpG context cytosines. This significantly reduces the
output size in most genomes. Note that using this option does not
merge data as symmetric CpGs.

```txt
-v, -verbose
```
Report more information while the program is running.
