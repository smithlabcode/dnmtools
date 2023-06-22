# roi - average methylation level in each region of interest

## Synopsis
```shell
$ dnmtools roi [OPTIONS] <intervals.bed> <input.meth>
```

## Description

One of the most common analysis tasks is to compute the average
methylation level through each of a set of genomic intervals. The
`roi` command does this. It reports weighted mean methylation in
each interval, along with read counts contributing to this mean. It
can also report an unweighted mean and a fraction of sites "called"
methylated in each interval. More details on these quantities can be
found in the documentation for the `levels` command.

The `roi` command requires two input files. The first is a
sorted [counts output file](../counts),
i.e. `input.meth` in the example above. This file provides data for
every site, either a cytosine or CpG, that is of interest. The second
input file (`intervals.bed`) specifies the genomic intervals in which
methylation statistics should be summarized. If either file is not
sorted by (chrom,end,start,strand) it can be sorted using the
following command:
```shell
$ LC_ALL=C sort -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o input-sorted.meth input.meth
```

The intervals must be specified as a BED format file, and these can be
sorted using [bedtools
sort](https://bedtools.readthedocs.io/en/latest/content/tools/sort.html).
Note that the `roi` command accepts intervals in two different BED
formats: (1) 6-column BED format, which may have more than 6 columns,
but requires the first 6 columns to match the specification, or (2)
3-column BED format.

From there, the `roi` command can be run as follows:
```shell
$ dnmtools roi -o output.bed regions.bed input-sorted.meth
```

The default output format is a 6-column BED format file, with the
"score" column (column 5 in the 6-column BED format) of each interval
now indicating the methylation level for the interval. If the input
intervals are 6-column BED format, then the default output only
differs in the "score" column. If the input intervals are a 3-column
BED file, then a name (auto-numbered and starting with "X") will be
assigned to each interval and all will be assumed on the positive
strand.

The value reported in the "score" column of the output file is, by
default, the weighted mean methylation through the interval. This is
calculated the same way as explained for the `levels` command.  Users
may specify two other values (the unweighted mean and the "fractional"
methylation) with a command line argument. If the input regions are
as follows:
```txt
chr1  3011124  3015902
chr1  3015904  3016852
chr1  3017204  3017572
chr1  3021791  3025633
chr1  3026050  3027589
```
Then the output might appear like:
```txt
chr1    3011124 3015902 X0  0.458647    +
chr1    3015904 3016852 X1  0.866667    +
chr1    3017204 3017572 X2  0.946429    +
chr1    3021791 3025633 X3  0.938038    +
chr1    3026050 3027589 X4  0.927007    +
```
If the input were given as 6-column BED format:
```txt
chr1    3011124 3015902 REGION_A    0   +
chr1    3015904 3016852 REGION_B    0   +
chr1    3017204 3017572 REGION_C    0   -
chr1    3021791 3025633 REGION_D    0   -
chr1    3026050 3027589 REGION_E    0   -
```
Then the output might be as follows:
```txt
chr1    3011124 3015902 REGION_A    0.458647    +
chr1    3015904 3016852 REGION_B    0.866667    +
chr1    3017204 3017572 REGION_C    0.946429    -
chr1    3021791 3025633 REGION_D    0.938038    -
chr1    3026050 3027589 REGION_E    0.927007    -
```

If additional information is desired in the methylation summary for
each interval, users can request this with a command line argument.
The additional information is all 3 varieties of "levels" (weighted,
unweighted and fractional), along with related integer counts through
the intervals. These are as follows, appearing after the first 6
columns of the 6-column BED format:

 * (7)  weighted mean methylation
 * (8)  unweighted mean methylation
 * (9)  fractional methylation
 * (10) number of CpGs in the region
 * (11) number of CpGs covered at least once
 * (12) number of observations in reads indicating methylation
 * (13) total number of observations from reads in the region

The weighted mean methylation level is then (12) divided by (13)
above. Example output might look like this:
```txt
chr1    3011124 3015902 REGION_A    0.458647    +   0.458647    0.448519    0.459302    172 172 915 1995
chr1    3015904 3016852 REGION_B    0.866667    +   0.866667    0.863492    1   6   6   39  45
chr1    3017204 3017572 REGION_C    0.946429    -   0.946429    0.954545    1   11  11  106 112
chr1    3021791 3025633 REGION_D    0.938038    -   0.938038    0.935424    1   109 109 1090    1162
chr1    3026050 3027589 REGION_E    0.927007    -   0.927007    0.918554    0.923077    13  13  127 137
```

Clearly if there are no reads mapping in a region, then the
methylation level will be undefined. By default these are indicated
with a value of "NA" in the output. This violates the 6-column BED
format convention of using a numerical score in the 5th column
(although strictly there are additional constraints on the possible
numerical values of the score column in BED format). So the option
`-N` for "numeric" will exclude such intervals in the output. It tends
to be rare that downstream analysis will not be robust to the "NA"
value, for example any processing in `R` can interpret the NA. Note:
when using the `-N` flag, if several methylation files are used with
the same set of genomic intervals, the output files may not have the
same number of intervals.

By default `roi` performs a binary search on the methylation file to
find sites (C or CpG) associated with each interval and loads only
those sites, which often saves both time and memory. However, it is
routine to compute the average methylation state across a large number
of target regions (e.g., fixed-size bins tiling the genome). When the
number of intervals is sufficiently large, it can be faster to load
the entire methylation file first. The `-L` option loads all lines of
the methcounts file into memory, which saves time at the expense of an
increased memory requirement.

## Options
```txt
-o, -output
```
The name of the output file. If no file name is provided, the output
will be written to standard output. Due to the size of this output, a
file should be specified unless the output will be piped to another
command or program. The output file contains genomic intervals in BED
format, with intervals corresponding to those provided as input.

```txt
-N, -numeric
```
print numeric values only (not NAs)
```txt
-L, -preload
```
Load all CpG sites
```txt
-s, -sort
```
sort data if needed
```txt
-l, -level
```
the level to report as score column in bed format output (w, u or f),
corresponding to weighted, unweighted or fractional methylation (default: w)

```txt
-M, -more-levels
```
report more methylation information
```txt
-v, -verbose
```
print more run info
