# fastlift - Mapping methylomes between species

## Synopsis
```shell
$ dnmtools fastlift -i <input.index> -f <input.from> -t <output.to>
```

Mapping methylomes between species builds on the
[liftOver tool](http://genome.ucsc.edu/cgi-bin/hgLiftOver) provided by
[UCSC Genome Browser](https://genome.ucsc.edu). However it is time
consuming to convert each methcounts output file from one assembly
to another using the UCSC liftOver tool, given that they all should
have the same locations but different read counts. Therefore, we use
liftOver to generate an index file between two assemblies, and provide
the `fast-liftover` tool.  Suppose we have downloaded the `liftOver` tool
and the chain file `mm9ToHg19.over.chain.gz` from the UCSC Genome
Browser website. If we have a methcounts file `mm9.meth` of
CpG sites or all cytosines in mm9.  Entries in  `mm9.meth`
look like

```txt
chr1  3005765  +  CpG  0.166667   6
chr1  3005846  +  CpG  0.5        10
chr1  3005927  +  CpG  0          9
```

We would like to lift it over to the human genome hg19, and generate
an index file `mm9-hg19.index` to facilitate later lift-over
operations from mm9 to hg19, and keep a record of unlifted mm9
cytosine positions in the file `mm9-hg19.unlifted`. First, convert the
[counts](../counts) file `mm9.meth` to the
BED file `mm9-cpg.bed` file for liftOver using the following command.

```shell
$ awk '{print $1"\t"$2"\t"$2+1"\t",$1":"$2":"$2+1":+\t0\t+"}' mm9.meth >mm9-cpg.bed
```

The output file `mm9-cpg.bed` should look like this:

```txt
chr1  3005765   3005766  chr1:3005765:3005766:+  0  +
chr1  3005846   3005847  chr1:3005846:3005847:+  0  +
chr1  3005927   3005928  chr1:3005927:3005928:+  0  +
```

Note that the fourth column is the genomic location data linked with
colons.

Then, run UCSC Genome Browser tool `liftOver` as follows:

```shell
$ liftOver mm9-cpg.bed mm9ToHg19.over.chain.gz mm9-hg19.index mm9-hg19.unlifted
```

The generated index file `mm9-hg19.index` will be a BED format file in
hg19 coordinates, with entries like

```txt
chr8    56539820        56539821  chr1:3005765:3005766:+        0       -
chr8    56539547        56539548  chr1:3005846:3005847:+        0       -
chr8    56539209        56539210  chr1:3005927:3005928:+        0       -
```

where the 4th column contains the genomic position of the cytosine
site in mm9 coordinates.

Next, convert the file `mm9-hg19.index` to a tab-separated input to be
passed onto the fast-liftover tool as follows.

```shell
$ tr ':' '\t' <mm9-hg19.index | awk '{print $4"\t"$5"\t"$1"\t"$2"\t"$9}' >mm9-hg19-fastliftover.index
```

After the index file is converted, we can use the `fast-liftover`
program on any mm9 methcounts file to lift it to hg19:

```shell
$ dnmtools fastlift -i mm9-hg19-fastliftover.index -f mm9.meth -t hg19-lift.meth
```

The `-p` option should be specified to report positions on the
positive strand of the target assembly.  Before using the lifted
methcounts file, make sure it is sorted properly.

```shell
$ LC_ALL=C sort -k1,1 -k2,2g -k3,3 hg19-lift.meth -o hg19-lift-sorted.meth
```

## Options
```txt
 -i, -indexfile
```
index file [required]
```txt
 -f, -from
```
Original file [required]
```txt
 -t, -to
```
 Output file liftovered [required]
```txt
 -u, -unlifted
```
(optional) File for unlifted sites
```txt
 -p, -plus-strand
```
 (optional) Report sites on + strand
```txt
 -v, -verbose
```
print more run info to STDERR as the program runs

