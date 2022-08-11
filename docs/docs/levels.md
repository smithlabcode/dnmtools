# levels - Methylation summary statistics

## Synopsis

```
$ dnmtools levels [OPTIONS] <input.meth>
```

## Description

The `levels` program computes statistics for the output of
[counts](../../analysis/counts).  Sample output is below. It
computes the total fraction of cytosines covered, the fraction of
cytosines that have mutated away from the reference, and coverage
statistics for both CpGs and all cytosines.

For CpG sites, coverage number reflects taking advantage of their
symmetric nature and merging the coverage on both strands. For CpG
coverage minus mutations, we remove the reads from CpG sites deemed
to be mutated away from the reference. It also computes average
methylation in three different ways, described in Schultz et al.
(2012). This program should provide flexibility to compare methylation
data with publications that calculate averages different ways and
illustrate the variability of the statistic depending on how it is
calculated. The sample output below only shows the results for
cytosines and CpGs in the sample, but similar statistics are provided
for symmetric CpGs and cytosines within the CHH, CCG, and CXG
contexts.

```
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

To run the levels program, execute
```
$ dnmtools levels -o output.levels input.meth
```

## Options

```
 -o, -output
```
output file in YAML format. (default: stdout)

```
 -a, -alpha
```
alpha for confidence interval  (default: 0.95)

```
 -v, -verbose
```
print more run info to STDERR while the program is running.
