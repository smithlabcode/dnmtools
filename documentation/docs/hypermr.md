# hypermr - Detecting hypermethylated regions

## Synopsis
```
$ dnmtools hypermr [OPTIONS] <input.meth>
```

## Description

The plant genomes, exemplified by *A. thaliana*, are devoid of DNA
methylation by default, with genic regions and transposons being
hyper-methylated, which we termed HyperMRs to stress their difference
from hypo-methylated regions in mammalian methylomes. DNA methylation
in plants has been associated with expression regulation and
transposon repression, and therefore characterizing HyperMRs is of
much biological relevance. In addition to plants, hydroxymethylation
tends to appear in a small fraction of the mammalian genome, and
therefore it makes sense to identify hyper-hydroxymethylated regions.

The first kind of HyperMR analysis involves finding continuous blocks
of hyper-methylated CpGs with the hmr program. Since the
[hmr](../hmr) program is designed to find hypo-methylated
regions, one can use it to identify HyperMRs by inverting the
methylation levels in the methcounts output file as follows:

```
$ awk '{$5=1-$5; print $0}' input.meth > input_inverted.meth
```

Next one may use the hmr program to find "valleys" in the inverted
Arabidopsis methylome, which are the hyper-methylated regions in the
original methylome. The command is invoked as below

```
$ dnmtools hmr -o output.hmr input_inverted.meth
```

This kind of HyperMR analysis produces continuous blocks of
hyper-methylated CpGs. However in some regions, intragenic regions in
particular, such continuous blocks of hyper-methylated CpGs are
separated by a few unmethylated CpGs, which have distinct sequence
preference when compared to those CpGs in the majority of unmethylated
genome. The blocks of hyper-methylated CpGs and gap CpGs together form
composite HyperMRs. The hypermr program, which implements a
three-state HMM, is used to identify such HyperMRs. Suppose the
[counts](../counts) output file is Col0 Meth.bed, to
find HyperMRs from this dataset, run

```
$ dnmtools hypermr -o output.hypermr input.meth
```

The output file is a 6-column
[BED](https://en.wikipedia.org/wiki/BED_(file_format))  file. The
first three columns give the chromosome, starting position and ending
position of that HyperMR.  The fourth column starts with the `hyper:`,
followed by the number of CpGs within this HyperMR. The fifth column
is the accumulative methylation level of all CpGs. The last column
indicates the strand, which is always +.

Lastly, it is worth noting that plants exhibit significantly more
methylation in the non-CpG context, and therefore inclusion of non-CpG
methylation in the calling of hyper-methylated regions could possibly
be informative. We suggest separating each cytosine context from the
methcounts output file as illustrated in the previous section (via
grep) and calling HyperMRs separately for each context.

## Options

```
 -o, -out
```
output BED file (default: STDOUT)

```
 -s, -scores
```
output file for posterior scores

```
 -t, -tolerance
```
tolerance (default: 0)

```
 -d, -desert
```
maximum distance between covered CpGs in HyperMR (default: 1000)

```
 -i, -itr
```
max number of iterations (default: 100)

```
 -V, -viterbi
```
Use Viterbi decoding

```
 -M, -min-meth
```
min cumulative methylation level in HypeMR (default: 4)
```
 -v, -verbose
```
print more run info to STDERR while the program is running.
```
 -P, -params-in
```
HMM parameters input file
```
 -p, -params-out
```
HMM parameters output file

