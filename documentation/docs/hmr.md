# hmr - hypomethylated regions

## Synopsis
```shell
$ dnmtools hmr [OPTIONS] <input.meth>
```

## Description
The distribution of methylation levels at individual sites in a
methylome (either CpGs or non-CpG Cs) almost always has a bimodal
distribution with one peak low (very close to 0) and another peak high
(close to 1). In most mammalian cells, the majority of the genome has
high methylation, and regions of low methylation are typically more
interesting.  These are called hypo-methylated regions (HMRs). In
plants, most of the genome has low methylation, and it is the high
parts that are interesting. These are called hyper-methylated regions.

For stupid historical reasons in the Smith lab, we called both of
these kinds of regions HMRs. One of the most important analysis tasks
is identifying the HMRs, and we use the `hmr` program for this. The
`hmr` program uses a hidden Markov model (HMM) approach using a
Beta-Binomial distribution to describe methylation levels at
individual sites while accounting for the number of reads informing
those levels. `hmr` automatically learns the average methylation
levels inside and outside the HMRs, and also the average size of those
HMRs.

## Requirements on the data

We typically like to have about 10x coverage to feel very confident in
the HMRs called in mammalian genomes, but the method will work with
lower coverage. Coverage can be calculated using the
[levels](../levels) program, and is summarized in the
`mean_depth_covered` statistic under `cpg_symmetric` group.

If reads have low coverage, the boundaries of HMRs will be less
accurate , but overall most of the HMRs will probably
be in the right places if you have coverage of 5-8x (depending on the
methylome). Boundaries of these regions are totally ignored by
analysis methods based on smoothing or using fixed-width windows.

## Typical mammalian methylomes

Running `hmr` requires a file of methylation levels formatted like the
output of the [counts](../counts). For calling HMRs in
mammalian methylomes, we suggest only considering the methylation
level at CpG sites, as the level of non-CpG methylation is not usually
more than a few percent. The required information can be extracted and
processed by using [sym](../sym).

## Output

The output will be in
[BED](https://en.wikipedia.org/wiki/BED_(file_format)) format, and the
indicated strand (always positive) is not informative. The name column
in the output will just assign a unique name to each HMR, and the
score column indicates how many CpGs exist inside the HMR. Each time
the `hmr` is run it requires parameters for the HMM to use in
identifying the HMRs. We usually train these HMM parameters on the
data being analyzed, since the parameters depend on the average
methylation level and variance of methylation level; the variance
observed can also depend on the coverage. However, in some cases it
might be desirable to use the parameters trained on one data set to
find HMRs in another. The option `-p` indicates a file in which the
trained parameters are written, and the argument `-P` indicates a file
containing parameters (as produced with the `-p` option on a previous
run) to use:

```shell
$ dnmtools hmr -p params.txt -o output.hmr input.meth
```

In the above example, the parameters were trained on the ESC
methylome, stored in the file human `params.txt` and then used to
find HMRs in the input methylome. This is useful if a particular
methylome seems to have very strange methylation levels through much
of the genome, and the HMRs would be more comparable with those from
some other methylome if the model were not trained on that strange
methylome.

### Partially methylated regions (PMRs)

The `hmr` program also has the option of directly identifying partially
methylated regions (PMRs), not to be confused with [partially
methylated domains](../pmd).  These are contiguous intervals
where the methylation level at individual sites is close to 0.5.  This
should also not be confused with regions that have allele-specific
methylation (ASM) or regions with alternating high and low methylation
levels at nearby sites.  Regions with ASM are almost always among the
PMRs, but most PMRs are not regions of ASM. The hmr program is run
with the same input but a different optional argument to find PMRs:

```shell
$ dnmtools hmr -partial -o human_esc.pmr human_esc.meth
```

## Converting HMR files to UCSC genome browser tracks

You might want to create bigBed browser tracks for HMRs.  The same
procedure also works for [AMRs](../amrfinder),
[PMDs](../pmd), or [DMRs](../dmr). To do so, follow these
steps:

 * (1) Download the bedToBigBed program from the UCSC Genome Browser
   [directory of binary utilities](http://hgdownload.cse.ucsc.edu/admin/exe/).
 * (2) Use the fetchChromSizes script from the same directory to
   create the \*.chrom.sizes file for the UCSC database you are working
   with (e.g. hg19). Note that this is the file that is referred to as
   `hg19.chrom.sizes` in step 3.
 * (3) Modify and use the following commands: PMDs, HMRs and AMRs may
   have a score greater than 1000 in the 5th column, in which case
  `bedToBigBed` will output an error. Also,  HMR file `input.bed` may have
   non-integer score in their 5th column.  The following script rounds
   the 5th column and prints 1000 if the score is bigger than 1000:
```shell
$ awk -v OFS="\t" '{if ($5>1000) print $1,$2,$3,$4,"1000"; else print $1,$2,$3,$4,int($5) }' input.bed > input.tobigbed
```
In the above command, since the HMRs are not stranded, we do not print
the 6th column. Keeping the 6th column would make all the HMRs appear
as though they have a direction – but it would all be the + strand. To
maintain the 6th column, just slightly modify the above awk command:
```shell
$ awk -v OFS="\t" '{if($5>1000) print $1,$2,$3,$4,"1000",$6; else print $1,$2,$3,$4,int($5),$6 }' human_esc.hmr > human_esc.hmr.tobigbed
```
 * (4) Generate the .bb track using the command below:
```shell
$ bedToBigBed input.tobigbed hg19.chrom.sizes output.bb
```

## Options

```txt
 -o, -out
```
output file (default: STDOUT)
```txt
 -d, -desert
```
maximum distance between covered CpGs in HMR (default: 1000)

```txt
 -i, -itr
```
max number of iterations (default: 100)
```txt
 -v, -verbose
```
print more run info to STDERR while the program is running.
```txt
 -partial
```
identify PMRs instead of HMRs
```txt
 -post-hypo
```
output file for single-CpG posterior hypomethylation probability (default: none)

```txt
 -post-meth
```
output file for single-CpG posteiror methylation probability (default: none)

```txt
 -P, -params-in
```
HMM parameter file (override training step)
```txt
 -p, -params-out
```
write HMM parameters to this file (default: none)
```txt
 -s, -seed
```
specify random seed (default: 408)
