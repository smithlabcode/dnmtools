# sym - collapse counts for symmetric CpGs sites

## Synopsis
```console
$ dnmtools sym [OPTIONS] <input.meth>
```

## Description

Many of our tools were designed for data vertebrate species. In these
species, the methylation levels at CpG sites tends to be symmetric,
the same on each strand. Of course there are exceptions. But in many
analysis settings, combining data from both strands for the same CpG
site is a good idea. Assume you have output from
[counts](../counts). The `sym` command will merge data on both strands
for each CpG site. It takes files having the same format as output by
`counts` with either all cytosines or CpGs only (generated with `-n`
option when running `counts`).
```console
$ dnmtools sym -o human_esc_CpG.meth human_esc.meth
```
The above command will merge all CpG pairs while also discarding sites
with an indication that the CpG has mutated. Note that as long as one
site of the pair is mutated, the pair is discarded. This is the
default mode. If you want to keep those mutated sites, run the
following:
```console
$ dnmtools sym -m -o human_esc_CpG.meth human_esc.meth
```

## Options

```txt
-o, -output
```
The name of the output file (default: stdout). The format is
the same as output by [counts](../counts).

```txt
-m, -muts
```
Include mutated CpG sites among the output, i.e. entries with an "x"
terminating the fourth column of each line of input.

```txt
-v, -verbose
```
Report more information while the program is running.
