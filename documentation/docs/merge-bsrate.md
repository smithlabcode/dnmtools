# merge-bsrate - Combine bisulfite conversion rate statistics files

## Synopsis
```shell
$ dnmtools merge-bsrate [OPTIONS] <file-1.bsrate> <file-2.bsrate> ...
```

## Description

Given several bisulfite conversion summary statistics generated using
the [bsrate](../bsrate) program, the `merge-bsrate` utility
combines their information. This is usually useful if your dataset has
been split into multipe files and processed in parallel, after which
one would like to combine the summaries of separate runs.

## Options

```txt
 -o -output
```
output file (default : STDOUT)

```txt
 -v -verbose
```
print more run info to STDERR as the program runs.
