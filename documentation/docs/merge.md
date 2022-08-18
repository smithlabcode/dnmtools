# merge - Combine counts files

## Synopsis
```console
$ dnmtools merge [OPTIONS] <file1.meth> <file2.meth> ...
```

## Description

This important command does two things, both of which can be called
"merging" of methylation files. The two behaviors are in one command
because the internal work is the same, but the output and the
motivation differs quite a bit between the two uses of merge.

(1) The merge command can take a set of [counts](../counts) output
files and combine them into one file. There are several reasons to do
this. One example is when technical replicates are performed, and
initially analyzed separately. But later these need to be combined for
downstream analysis. If the counts output files are already available,
then merging them using the merge command is easier than re-doing the
analysis with all replicates together. Similarly, if some data is
generated earlier, and more produced later from the same sample,
merging in this way avoids repeating some of the early analysis
effort.

Suppose you have the three counts output files from
three different replicates: `R1.meth`, `R2.meth` and
`R3.meth`. To merge those individual methcounts files, execute:
```console
$ dnmtools merge -o combined.meth R1.meth R2.meth R3.meth
```
The command can handle an arbitrary number of files, and the files
do not need to have the same number of lines/sites. The merge command
does assume that the sorted order of chromosomes within each input
file is consistent, and for each input file, within a chromosome
all sites appear in increasing order.

(2) The merge command can take a set of counts output files and
combine them as a table that contains all the same information. The
table format is helpful if subsequent analyses are to be done using a
data table, for example a data frame in R. When producing this tabular
format, merge allows the user to select whether the desired output is
in counts (both the count of methylated and unmethylated reads for
every site) or as the fractions. If the fraction would have involved
division by zero, then "NA" is written in the output. But this
behavior can be controlled with the command line options.

Suppose you have 4 different methylomes, two replicates from wild
type and two from a mutant. You want to create a table of
methylation information so you can analyze the data in R as a data
frame. The merge command can help by pasting the data as a table,
ensuring a consistent order and filling in any missing values from
the individual input files:
```console
$ dnmtools merge -o table.txt wt1.meth wt2.meth mut1.meth mut2.meth
```

The file `table.txt` is not in the same format as the input files,
since those each have exactly 6 columns. The output has one column as
the row names. Then it will have two columns for each of the input
files, one with the count of total reads, and one with the count of
reads indicating methylation. Here is what the output might look like:
```txt
                wt1_R   wt1_M   wt2_R   wt2_M   mut1_R   mut1_M mut2_R  mut2_M
chr1:108:+:CpG  9       6       10      8       2       2       2       1
chr1:114:+:CpG  17      7       10      0       5       1       9       1
chr1:160:+:CpG  12      8       10      5       15      14      13      6
chr1:309:+:CpG  1       1       1       0       12      8       2       1
```

Note: Currently the output from merge may not be immediately
compatible with [radmeth](../radmeth), since the column headings might
not in the expected format. An option for this has been added in the
most recent code.

## Options

```txt
-o, -output
```
output file as [counts](../counts) format (default: stdout)

```txt
-h, -header
```
Print a header given by the input string at the top of the file
(ignored for tabular)

```txt
-t, -tabular
```
Output is in table format.

```txt
-remove
```
Suffix to remove from filenames when making column names for tabular
format. If not specified, suffix including from final dot is removed.

```txt
-s, -suff
```
Column name suffixes, one for total reads and one for methylated
reads, to be separated from sample name with underscore in the header
for tabular format output.

```txt
 -f, -fractional
```
Output table will give fractions (requires `-tabular`).

```txt
 -r, -reads
```
Minimum number of reads required when using the `-f` flag (default: 1)

```txt
-ignore
```
Ignore sorting. Do not attempt to determine chromosome
order. Lexicographic order on chromosome names will be assumed.

```txt
 -v, -verbose
```
Report more information while the program is running.
