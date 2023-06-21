# Contexts for cytosine methylation

We only consider DNA methylation at cytosines, and the "contexts"
handled in dnmtools include the following:

* `CpG` sites: These are dinucleotides, so must span 2 positions. Each
  has two cytosines, one on each strand. Certain tools will work with
  "symmetric" CpG sites, which are associated with data from cytosines
  on both strands. This is because, if we exclude mutations and
  certain specific cases, the methylation status of the two cytosines
  at a CpG site tends to be the same. It is wise to check this
  assumption before depending on it too heavily. In a file listing
  locations of all CpG sites, including both strands, then CpG sites
  should appear in pairs at consecutive positions on the chromosome
  and with opposite strands. The "p" in the CpG symbol provides no
  additional information if we have strand and always assume 5'->3',
  but we use it for historical reasons.

* `CHH` sites: These are cytosines followed by two sites, neither
  being a guanine (H indicates "not G"). These cannot be symmetric,
  and have unique methylation properties in some plants.

* `CHG` sites: These are cytosines followed a non-guanine, and then by
  a guanine. These sites can considered "symmetric", with respect to
  their C and G, but none of our analysis tools currently consider
  them as such. Certain methylation properties at CHG sites, in some
  plants for example, have been found to differ from those having CpG
  and CHH contexts.

* `CCG` sites: These are a subset of the CHG sites, but they also
  overlap CpG sites (they include a CpG site) and there has been
  suggestion that their properties may be unique. We therefore mark
  these sites in some of our analysis tools. Although these could be
  considered symmetric, the complementary site would be a CGG, which
  has only one cytosine.

* `CXG` sites: These are also a subset of the CHG sites. We invented
  this notation, and it was a poor choice. A better choice would have
  been CWG, as W is IUPAC for A or T. The CAG and CTG trinucleotides
  are common in mammals, and likely are related to mutation from CCG
  sites on one strand or the other. In a future version of our
  software we hope to change this symbol.

Depending on the command within dnmtools, the above "contexts" may be
used to form exclusive categories or just as subsets of sites.

Suppose we have the following DNA sequence as reference, with
positions numbered below the sequence and a `+` or `-` symbol to
indicate strand for each cytosine:
```txt
ATGCCCGCAAGGTCTG
  -+++-+  -- + -
0123456789012345
0000000000111111
```
The "context" for each cytosine, along with the position
and strand, would be as follows:
```txt
 2   -   CHH
 3   +   CHH
 4   +   CCG
 5   +   CpG
 6   -   CpG
 7   +   CHH
10   -   CHH
11   -   CHH
13   +   CXG
15   -   CXG
```
If we were not distinguishing the `CCG` from the `CXG`, we the
contexts would be:
```txt
 2   -   CHH
 3   +   CHH
 4   +   CHG
 5   +   CpG
 6   -   CpG
 7   +   CHH
10   -   CHH
11   -   CHH
13   +   CHG
15   -   CHG
```

The following table should help explain which triples of nucleotides
are counted towards each context. Each of the triples begins with a C,
and in our formats, this is the cytosine where the methylation level
or state is in question.

| Trip | CpG | CHH | CHG | CWG | CCG |   |
|------|-----|-----|-----|-----|-----|---|
| CAA  |     | 1   |     |     |     |   |
| CAC  |     | 1   |     |     |     |   |
| CAG  |     |     | 1   | 1   |     | * |
| CAT  |     | 1   |     |     |     |   |
| CCA  |     | 1   |     |     |     |   |
| CCC  |     | 1   |     |     |     |   |
| CCG  |     |     | 1   |     | 1   | * |
| CCT  |     | 1   |     |     |     |   |
| CGA  | 1   |     |     |     |     |   |
| CGC  | 1   |     |     |     |     |   |
| CGG  | 1   |     |     |     |     |   |
| CGT  | 1   |     |     |     |     |   |
| CTA  |     | 1   |     |     |     |   |
| CTC  |     | 1   |     |     |     |   |
| CTG  |     |     | 1   | 1   |     | * |
| CTT  |     | 1   |     |     |     |   |

The traditional contexts considered are the CpG, the CHH and the
CHG. The CHH and CHG have been of more interest in plants, especially
Arabidopsis. Together the CpG, CHH and CHG contexts cover all
trinucleotides that start with C, and partition the trinucs
unambiguously. These are the first 3 columns above. The CWG and CCG
are of interest mostly because the CWG is so important in vertebrate
species (and more-so for mammalia) where a combination of
deamination-induced loss of the middle cytosine of a CCG and tandem
expansion of CAG/CTG repeats (including within human populations) has
led to a relative abundance of the CWG in important places in the
genome. The CWG may also be called "symmetric" in the same way as a
CHG. However, in a strict sense, if one calls every CHG symmetric,
then it might include CCG on one strand, with CGG on the other, and
the CGG would not be a CHG.

The above table does not mention the CXG; we plan to remove the CXG
and only include contexts from the above table.
