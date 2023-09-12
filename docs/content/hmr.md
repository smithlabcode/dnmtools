# hmr - hypomethylated regions

## Synopsis
```console
$ dnmtools hmr [OPTIONS] <input.meth>
```

## Description
This command identifies "hypomethylated regions" which we abbreviate
as HMRs. The valleys are identified using a 2-state hidden Markov
model with a beta-binomial emission distribution for each CpG
site. This distribution accounts for methylation and random changes to
coverage along the genome. `hmr` automatically learns the average
methylation levels inside and outside the HMRs, and also the average
size of those HMRs.

The input to `hmr` is a file with the output format from the `counts`
command. The output is a BED format file containing the HMRs as
genomic intervals.

For mammalian healthy somatic primary cells, these are the most
important functional features in the methylome. They tend to
correspond to promoters and enhancers (and other regions of possible
regulatory activity) that are functional/active, are accessible and
poised for function, or that were so in a progenitor cell. Although
the latter two categories are important, in most somatic primary cells
the HMRs mark active regulatory regions. From a global view of the
methylome, these are the valleys in an otherwise high background
methylation level. In a typical sample of healthy somatic primary
cells, you should expect to find between roughly 40k-100k HMRs, and
their mean size should be 1.5 kbp to 3 kbp. If your results deviate
too much from this, then you should consider whether it makes sense to
identify HMRs in your sample (e.g., if it's a cancer sample, or
immortalized, then the HMRs will be obscured by other features). I
have never observed HMRs defined in this way in species outside
vertebrates. For example, Arabidopsis has a low background methylation
level punctuated by peaks, rather than valleys (and the difference
isn't as simple as subtracting the methylation level from 1.0).

## Requirements on the data

Running `hmr` requires a file of methylation levels formatted like the
output of the [counts](../counts). For identifying HMRs in mammalian
methylomes, use the symmetric CpG methylation levels. This is obtained
by using the [sym](../sym) command after having used the
[counts](../counts) command.

We typically like to have about 10x coverage to feel very confident in
the HMRs called in mammalian genomes, but the method will work with
lower coverage. Coverage can be calculated using the
[levels](../levels) command, and is summarized in the
`mean_depth_covered` statistic under `cpg_symmetric` group.

If reads have low coverage, the boundaries of HMRs will be less
accurate, but overall most of the HMRs will probably be in the right
places if you have coverage of 5-8x (depending on the methylome).
Boundaries of these regions are ignored by analysis methods based on
smoothing or using fixed-width windows, so you will get better
precision on boundaries (and accuracy overall) using `hmr`.

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
```console
$ dnmtools hmr -p params.txt -o output.hmr input.meth
```
Above the output file has the extension `.hmr` but this doesn't
matter. The format of the output is 6-column BED.

In the above example the trained parameters are stored in the file
`params.txt` but are also used to find HMRs in the input
methylome. Storing these parameters can be useful if a particular
methylome seems to have very strange methylation levels through much
of the genome, and the HMRs would be more comparable with those from
some other methylome if the model were not trained on that strange
methylome.

**Partially methylated regions (PMRs)**

The `hmr` command also has the option of directly identifying partially
methylated regions (PMRs), not to be confused with [partially
methylated domains](../pmd).  These are contiguous intervals
where the methylation level at individual sites is close to 0.5.  This
should also not be confused with regions that have allele-specific
methylation (ASM) or regions with alternating high and low methylation
levels at nearby sites.  Regions with ASM are almost always among the
PMRs, but most PMRs are not regions of ASM. The hmr command is run
with the same input but a different optional argument to find PMRs:
```console
$ dnmtools hmr -partial -o output.pmr input.meth
```

## Options

```txt
 -o, -out
```
The name of the output file. If no file name is provided, the output
will be written to standard output. Due to the size of this output, a
file name should be specified unless the output will be piped to
another command or program. The output file contains genomic intervals
in BED format.

```txt
 -d, -desert
```
The maximum distance between covered CpGs in HMR (default: 1000
bp). Beyond this distance, adjacent CpG sites will be considered part
of distinct HMRs, regardless of their methylation status.

```txt
 -i, -itr
```
The maximum number of iterations for learning parameters (default:
10).

```txt
 -v, -verbose
```
Report more information while the program is running.

```txt
 -partial
```
Identify PMRs instead of HMRs.

```txt
 -post-hypo
```
Output file for single-CpG posterior hypomethylation probability
(default: none). By default this is information is not reported, and
only reported if a file is specified here.

```txt
 -post-meth
```
Output file for single-CpG posteiror methylation probability (default:
none). By default this is information is not reported, and only
reported if a file is specified here.

```txt
 -P, -params-in
```
File containing existing parameters to use in the model (skip the
training step). This should be a file produced previously by the
hmr command using the `-p` parameter.

```txt
 -p, -params-out
```
File in which to write parameters learned during the current run.

```txt
 -s, -seed
```
A random number seed. Randomization is used in a shuffling step prior
to filering candidate HMRs. This parameter is typically only used for
testing (default: 408).

```txt
 -summary
```
Write the analysis summary to this file. The summary is not
reported unless a file is specified here.
