# bsrate - estimate bisulfite conversion rate

## Synopsis
```shell
$ dnmtools bsrate [OPTIONS] -c <chroms> <input.sam>
```

## Description

Unmethylated cytosines in DNA fragments are converted to uracil by
sodium bisulfite treatment. As these fragments are amplified, the
uracils are converted to thymines and so unmethylated Cs are
ultimately read as Ts (assuming no error). Despite its high fidelity,
bisulfite conversion of C to T does have some inherent failure rate,
depending on the bisulfite kit used, reagent concentration, time of
treatment, etc. These factors impact the success rate of the reaction.
Therefore, the bisulfite conversion rate, defined as the rate at which
unmethylated cytosines in the sample appear as Ts in the sequenced
reads, should be measured and should be very high (e.g. > 0.99) for
the experiment to be considered a success.

Measuring bisulfite conversion rate this way requires some kind of
control set of genomic cytosines not believed to be methylated. Three
options are (1) to spike in some DNA known not to be methylated, such
as a Lambda virus, (2) to use the human mitochondrial genome, which is
known to be unmethylated almost everywhere in nearly every cell type
(or chloroplast genomes), or (3) to use non-CpG cytosines which are
believed to be almost completely unmethylated in most mammalian
cells. In general the procedure is to identify the positions in reads
that correspond to these presumed unmethylated cytosines, then compute
the ratio of T over (C + T) at these positions. If the bisulfite
reaction is perfect, then this ratio should be very close to 1, and if
there is no bisulfite treatment, then this ratio should be close to 0.

The `bsrate` command will estimate the bisulfite conversion rate in
this way. Assuming method (3) from the above paragraph of measuring
conversion rate at non-CpG cytosines in a mammalian methylome, the
following command will estimate the conversion rate.

```shell
$ dnmtools bsrate -c /path/to/genome.fa -o output.bsrate input-sorted.sam
```

Note that we often use the output of [uniq](../uniq) to
reduce any bias introduced by differential PCR amplification as a
function of conversion. The `bsrate` command requires that the input
be sorted so that reads mapping to the same chromosome are consecutive
within the file and the input should have duplicate reads already
removed. The first several lines of the output might look like the
following:

```txt
OVERALL CONVERSION RATE = 0.994141
POS CONVERSION RATE = 0.994166  832349
NEG CONVERSION RATE = 0.994116  825919
BASE PTOT  PCONV PRATE   NTOT  NCONV NRATE   BTHTOT BTHCONV BTHRATE ERR ALL    ERRRATE
1    8964  8813  0.9831  9024  8865  0.9823  17988  17678   0.9827  95  18083  0.0052
2    7394  7305  0.9879  7263  7183  0.9889  14657  14488   0.9884  100 14757  0.0067
3    8530  8442  0.9896  8323  8232  0.9890  16853  16674   0.9893  98  16951  0.0057
4    8884  8814  0.9921  8737  8664  0.9916  17621  17478   0.9918  76  17697  0.0042
5    8658  8596  0.9928  8872  8809  0.9929  17530  17405   0.9928  70  17600  0.0039
6    9280  9218  0.9933  9225  9177  0.9948  18505  18395   0.9940  59  18564  0.0031
7    9165  9117  0.9947  9043  8981  0.9931  18208  18098   0.9939  69  18277  0.0037
8    9323  9268  0.9941  9370  9314  0.9940  18693  18582   0.9940  55  18748  0.0029
9    9280  9228  0.9944  9192  9154  0.9958  18472  18382   0.9951  52  18524  0.0028
10   9193  9143  0.9945  9039  8979  0.9933  18232  18122   0.9939  66  18298  0.0036
```

The above example is based on a very small number of mapped reads in
order to make the output fit here. The first thing to notice is that
the conversion rate is computed separately for each strand. The
information is presented separately because this is often a good way
to see when some problem has occurred in the context of paired-end
reads. If the conversion rate looks significantly different between
the two strands, then we would go back and look for a mistake that has
been made at an earlier stage in the pipeline. The first 3 lines in
the output indicate the overall conversion rate, the conversion rate
for positive strand mappers, and the conversion rate for negative
strand mappers. The total number of nucleotides used (e.g. all C+T
mapping over genomic non-CpG C’s for method (3)) is given for positive
and negative strand conversion rate computation, and if everything has
worked up to this point these two numbers should be very similar. The
4th line gives column labels for a table showing conversion rate at
each position in the reads. The labels PTOT, PCONV and PRATE give the
total nucleotides used, the number converted, and the ratio of those
two, for the positive-strand mappers. The corresponding numbers are
also given for negative strand mappers (NTOT, NCONV, NRATE) and
combined (BTH). The sequencing error rate is also shown for each
position, though this is an underestimate because we assume at these
genomic sites any read with either a C or a T contains no error.

The bisulfite conversion rate reported in MethBase is computed from
all non-CpG cytosines, and assumes zero non-CpG methylation. Because
this is likely not 100% true, it is a conservative estimate of
conversion rate. A better method for computing conversion rate is to
use an unmethylated spike-in or reads mapping to mitochondria, which
has been shown to be entirely unmethylated in most human tissues. To
use the mitochondria, you can extract mitochondrial reads from a SAM
file using `awk`:

```shell
$ awk '$1 /~/ ^@ || $3 == "chrM"' input.sam >input-chrM-only.sam
```

After completing bisulfite conversion rate analysis, remember to
remove any control reads not naturally occurring in the sample (lambda
virus, mitochondrial DNA from another organism, etc.) before
continuing your analysis. The output from two different runs of
`bsrate` can be merged using the command
[merge-bsrate](utils/merge-bsrate).

## Options

```txt
 -o, -output
```

The name of the output file (default: STDOUT).

```txt
 -c, -chrom
```

File or directory of chromosome sequences (FASTA format; .fa suffix)
[required]

```txt
 -N, -all
```

Count all Cs (including CpGs) when estimating bisulfite
conversion. This will only work if the estimate is made from sequences
known to be unmethylated, like with a spike-in.

```txt
 -seq
```

Use only reads that map to this chromosome (e.g. chrM).

```txt
 -v, -verbose
```

Print more information while the program is running.
