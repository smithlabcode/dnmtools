# Visualizing methylome data

Here we explain how to visualize data using the UCSC Genome
Browser. When we refer to the genome browser below, we mean the UCSC
kind.

## Single-site methylation levels

Here we are concerned with individual sites. These need not be CpG
sites -- the could be any/all cytosines, but we will assume they are
CpGs through our explanation.

To view the methylation level at individual CpG sites in a genome
browser, the data should be converted into bigWig format. The starting
point should be a "counts" file, as output from the
[counts](../counts) command. The bigWig format is intended for the
"wiggle" tracks, which shows information associated with individual
genomic positions, but in the bigWig format this information is
encoded concisely and is not for direct human viewing. The same
approach is used to build files that show the coverage at individual
CpG sites.

To create methylation level tracks or read coverage tracks, follow
these steps:

* Download the `wigToBigWig` program from the UCSC Genome Browser
  directory of [binaries](http://hgdownload.cse.ucsc.edu/admin/exe/).

* Use the `fetchChromSizes` script, from the same directory, to get
  the `.chrom.sizes` file for the database (reference genome) you are
  working with (e.g., hg38). Note: this is the file mentioned below as
  `hg19.chrom.sizes` for the hg19 reference genome.

* To create a bigWig format track for methylation levels at CpG
  sites, convert the symmetric methylation file ([counts](../counts)
  format) as follows:
```console
$ awk '{print $1,$2,$2+1,$5}' sample.meth | wigToBigWig /dev/stdin hg19.chrom.sizes sample.meth.bw
```
  In the command above, the first part selects the appropriate columns
  to generate bedgraph format, and then the second part converts this
  directly into a bigWig format file, which is not human-readable.

* To create a bigWig format track for read coverage at CpG sites, use the
  following command, which is very similar to the previous one above:
```console
$ awk '{print $1,$2,$2+1,$6}' sample.meth | wigToBigWig /dev/stdin hg19.chrom.sizes sample.reads.bw
```

If the `wigToBigWig` or `fetchChromSizes` programs are not
executable when downloaded, try the following:
```console
$ chmod +x wigToBigWig
$ chmod +x fetchChromSizes
```

## The identified features

This refers to the HMRs, the AMRs, the PMDs, and possibly the
HyperMRs. These are contiguous genomic intervals. It happens that for
an individual set of these features, as identified using dnmtools, no
two features will overlap. This fact isn't relevant here, though.

We will assume you want to make browser tracks for HMRs. The same
procedure also works for [AMRs](../amrfinder), [PMDs](../pmd), or
[DMRs](../dmr). To do so, follow these steps:

* Download the `bedToBigBed` program from the UCSC Genome Browser
  directory of [binaries](http://hgdownload.cse.ucsc.edu/admin/exe/).

* Use the `fetchChromSizes` script, from the same directory, to get
  the `.chrom.sizes` file for the database (reference genome) you are
  working with (e.g., hg38). Note: this is the file mentioned below as
  `hg19.chrom.sizes` for the hg19 reference genome.

* Modify and use the following commands: PMDs, HMRs and AMRs may have
  a score greater than 1000 in the 5th column, in which case
  `bedToBigBed` will output an error. Also, HMR file `sample.bed` may
  have a non-integer score in the 5th column. The following script
  rounds the 5th column and prints 1000 if the score is greater than
  1000:
```console
$ awk -v OFS="\t" '{if ($5>1000) print $1,$2,$3,$4,"1000"; \
                    else print $1,$2,$3,$4,int($5)}' sample.bed > sample.for_bigbed
```
  In the above command, since the HMRs are not stranded, we do not keep
  the 6th column. Keeping the 6th column would make all the HMRs appear
  as though they have a direction and they would all appear to be on the +
  strand. This would be visually misleading (and somewhat annoying).

* Generate the `.bb` track file using the command below:
```console
$ bedToBigBed sample.for_bigbed hg19.chrom.sizes output.bb
```
