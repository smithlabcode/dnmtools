# levels - global methylation summary statistics

## Synopsis

```console
$ dnmtools levels [OPTIONS] <input.meth>
```

## Description

The `levels` command computes global summary statistics for the output
of [counts](../counts). Example output is below. It computes multiple
summary statistics related to the quantity of data (e.g., coverage of
sites) and methylation (e.g., global average methylation). These
summary statistics are also provided by context. The contexts are not
exclusive categories, and include:

* cytosines, all of them, on either strand
* cpg sites, on either strand
* symmetric cpg sites (strands combined)
* the CHH context
* the CCG context
* the CXG context (we "invented" this one)

The summary statistics computed include:

* `total_sites` the total number of sites counted for this context
* `sites_covered` among the total above, those with at least one read
* `total_c` among the observations in reads, how many are C
* `total_t` among the observations in reads, how many are T
* `max_depth` the most coverage of any site for this context
* `mutations` number of sites for this context marked as mutated
* `called_meth` number of sites "called" methylated
* `called_unmeth` number of sites "called" unmethylated
* `mean_agg` the sum of methylation levels for all sites
* `coverage` total data informing on sites for this context
* `sites_covered_fraction` fraction of sites covered
* `mean_depth` among all sites, the mean coverage by reads
* `mean_depth_covered` among all covered sites, the mean coverage
* `mean_meth` the mean of the methylation levels for covered sites
* `mean_meth_weighted` the mean weighted by coverage
* `fractional_meth` the fraction of "called" sites "called" methylated

(If you want more information on these, please ask.)

Among the above, many are included because they are needed for
calculating the the "derived" statistics. For example, the `mean_agg`
is used in the denominator for `mean_meth`, where the denominator is
the number of covered sites. Why keep those raw statistics? Because
it's essential if two different `levels` output files are combined.

The final three values are the "levels" and are described in Schultz
et al. (2012):
```txt
"Leveling" the playing field for analyses of single-base resolution DNA methylomes
Schultz, Schmitz & Ecker (TIG 2012)
```

Note: the `fractional_meth` level we calculate is inspired but
different from the paper. What we are do is use a binomial test to
determine significantly hyper/hypomethylated sites, and only use the
subset of significant sites to calculate `fractional_meth` level.

This command should provide flexibility to compare methylation data
with publications that calculate averages different ways. The sample
output below only shows the results for cytosines and CpGs in the
sample, but similar output is generated for symmetric CpGs and
cytosines in the CHH, CCG, and CXG contexts.

```yaml
cytosines:
  total_sites: 1200559022
  sites_covered: 797100353
  total_c: 417377038
  total_t: 4048558428
  max_depth: 30662
  mutations: 3505469
  called_meth: 44229556
  called_unmeth: 750163257
  mean_agg: 4.40429e+07
  coverage: 4465935466
  sites_covered_fraction: 0.663941
  mean_depth: 3.71988
  mean_depth_covered: 5.60273
  mean_meth: 0.055254
  mean_meth_weighted: 0.093458
  fractional_meth: 0.055677
cpg:
  total_sites: 58803590
  sites_covered: 47880982
  total_c: 261807401
  total_t: 84403225
  max_depth: 30080
  mutations: 381675
  called_meth: 38861909
  called_unmeth: 7152004
  mean_agg: 3.69282e+07
  coverage: 346210626
  sites_covered_fraction: 0.814253
  mean_depth: 5.88758
  mean_depth_covered: 7.23065
  mean_meth: 0.771250
  mean_meth_weighted: 0.756208
```

You can run the `levels` command as follows:
```console
$ dnmtools levels -o output.levels input.meth
```

## Options

```console
-o, -output
```
Output file in YAML format (default: stdout).

```console
-a, -alpha
```
Alpha for confidence interval (default: 0.95).

```console
-v, -verbose
```
Report more information while the program is running.
