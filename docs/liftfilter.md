# liftfilter - merge lifted entries to the same position

## Synopsis
```shell
$ dnmtools liftfilter [OPTIONS] -o <output.meth> <input.meth>
```

## Description

The [fastlift](../fastlift) program may report multiple mm9
sites mapped to a same position in hg19.  In this situation, we may
either collapse read counts at those mm9 sites, or keep the data for
only one mm9 site. We can use the lift-filter program to achieve these
two options. Use

```shell
$ dnmtools liftfilter -o output-filtered.meth input.meth
```

to merge data from mm9 sites lifted to the same hg19 position. Use the
option `-u` to keep the first record of duplicated sites.

## Options

```txt
 -o, -output
```
Output processed methcount [required]
```txt
 -u, -unique
```
 keep unique sites
```txt
 -v, -verbose
```
print more information to STDERR as the program runs.
