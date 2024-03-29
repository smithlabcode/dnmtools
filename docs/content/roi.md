# roi - average methylation level in each region of interest

## Synopsis
```shell
$ dnmtools roi [OPTIONS] <intervals.bed> <input.counts>
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
i.e. `input.counts` in the example above. This file provides data for
every site, either a cytosine or CpG, that is of interest. The second
input file (`intervals.bed`) specifies the genomic intervals in which
methylation statistics should be summarized. If either file is not
sorted by (chrom,end,start,strand) it can be sorted using the
following command:
```shell
$ LC_ALL=C sort -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o input-sorted.counts input.counts
```
Note: As of v1.4.0, the sorted order of chromosomes/targets within these
files is not important, but the sites within each chromosome must
still be sorted.

The intervals must be specified as a BED format file, and these can be
sorted using [bedtools
sort](https://bedtools.readthedocs.io/en/latest/content/tools/sort.html).
Note that the `roi` command accepts intervals in two different BED
formats: (1) 6-column BED format, which may have more than 6 columns,
but requires the first 6 columns to match the specification, or (2)
3-column BED format.

*An important note about the input files:* several aspects of the
output for `roi` depend on the number of sites within each region of
interest. If the `.counts` file provided as input does not have all
the sites you might expect, for example if it is missing sites that
have been excluded from some earlier step in your pipeline, then the
results will be affected. We hope to make `roi` more robust to this
issue in the future, for example by accepting some information about
the reference genome to ensure that the numbers of sites are as
expected by the user.

From there, the `roi` command can be run as follows:
```shell
$ dnmtools roi -o output.bed regions.bed input-sorted.counts
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

**Performance**
There are two different methods for reading the `.counts` input file.
The default mode of assuming the counts file is sorted and doing
binary search on disk will only work if the file is not zipped. If the
file is not zipped, the user must specify `-L` to load all the counts
into memory. These two methods should produce identical output files,
but depending on the application one might be much faster than the
other. The numbers below are based on sampling random intervals from a
file with mean interval size of 2237 bp.  The time is measured in
seconds, and all work was done on local SSDs.

|     N | on disk | in memory |
|   ---:|     ---:|       ---:|
|    10 |  1.70 s |    9.55 s |
|   100 |  1.75 s |    9.43 s |
|  1000 |  2.05 s |    9.56 s |
| 10000 |  5.13 s |    9.64 s |
| 20000 |  8.49 s |    9.81 s |
| 30000 | 11.93 s |    9.83 s |
| 40000 | 15.20 s |   10.06 s |
| 50000 | 18.96 s |    9.95 s |
|    10 |  6.7 MB | 3079.7 MB |
| 50000 | 10.1 MB | 3084.0 MB |

The counts file for the above data had 30962770 lines, corresponding
to CpG sites in hg38 including those in the "extra" chromosomes.  Be
aware that if you work on a cluster or in the cloud, latency of the
disks might become a problem so the "in memory" option might be better
all the way. This depends on how your storage is configured. Even if
throughput is highly tuned, latency can cause major slowdown for the
"on disk" mode.

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
