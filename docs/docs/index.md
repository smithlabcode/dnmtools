# DNMTools

The `dnmtools` set of tools for analyzing DNA methylation from
high-throughput sequencing data, especially whole genome bisulfite
sequencing (WGBS), but also reduced representation bisulfite
sequencing (RRBS).

This set of tools is superseding our collection of tools in MethPipe.
The design of `dnmtools` differs, but it has all the same
functionality, with most tools in one-to-one correspondence.

At present, as we make the first release of `dnmtools`, the
documentation is structured very differently from that of MethPipe. We
will rebuild that documentation here, but it might take a few
iterations to stabilize. Most of the old documentation still applies
for tasks that do not directly involve our tools. If you are already a
user of our tools, please check the documentation here for the
individual commands within `dnmtools` to ensure you are aware of any
changes to command line arguments or file formats.

We now include a minimal pipeline (or set of pipelines) using
Snakemake and we expect this will expand.

# Assumptions

Our pipeline was designed to run in a cluster computing context, with
many processing nodes available, and a job submission system like PBS,
SGE or SLURM, but it is also possible to analyze methylomes from most
genomes (including human and mouse) in a local machine with at least
16 GB of RAM.  Typically the data we deal with amounts to a minimum of
100GB for a mammalian methylome at 10x coverage. Intermediate files
may cause this amount to more than double during execution of the
pipeline, and likely at the end of the pipeline the total size of
files will amount to almost double the size of the raw data.

It is critical that users are familiar with bisulfite sequencing
experiments, especially the bisulfite conversion reaction, and how
this affects what we observe in the sequenced reads.

Users are assumed to be somewhat familar with a command line.
Installing `dnmtools` might require working with some environment
variables. The data that `dnmtools` is designed for is really big, so
understanding your computer will be important.

# Citation information

If you find our programs helpful, please cite us. The majority of the
commands within `dnmtools` were introduced in their own
publications. However, several improvements and additions have been
made over the years. Below are the relevant citations.

**dnmtools**

Qiang Song, Benjamin Decato, Elizabeth E Hong, Meng Zhou, Fang Fang,
Jianghan Qu, Tyler Garvin, Michael Kessler, Jun Zhou, and Andrew D Smith.
*A reference methylome database and analysis pipeline to facilitate integrative and comparative epigenomics.*
PloS one, 8(12):e81148, 2013

**abismal**

Guilherme de Sena Brandine, and Andrew D Smith
*Fast and memory-efficient mapping of short bisulfite sequencing reads
using a two-letter alphabet.*
NAR Genomics and Bioinformatics, 3(4), lqab115, 2021.

**amrfinder**

Fang Fang, Emily Hodges, Antoine Molaro, Matthew Dean,
Gregory J Hannon, and Andrew D Smith.
*Genomic landscape of human allele-specific dna methylation.*
Proceedings of the National Academy of Sciences, 109(19):7332–7337, 2012

**mlml**

Jianghan Qu, Meng Zhou, Qiang Song, Elizabeth E Hong, and
Andrew D Smith.
*Mlml: Consistent simultaneous estimates of dna methylation and
hydroxymethylation.*
Bioinformatics, 29(20):2645–2646, 2013

**radmeth**

Egor Dolzhenko and Andrew D Smith.
*Using beta-binomiaregression for high-precision differen
tial methylation analysis in multifactor whole-genome bisulfite
sequencing experiments.*
BMC bioinformatics, 15(1):1–8, 2014

**pmd**

Benjamin E Decato, Jianghan Qu, Xiaojing Ji, Elvin Wagenblast,
Simon RV Knott, Gregory J Hannon, and Andrew D Smith.
*Characterization of universal features of partially methylated domains across tissues and
species.*
Epigenetics & Chromatin, 2020

Contacts and bug reports
========================

Andrew D. Smith
andrewds@usc.edu

Guilherme de Sena Brandine
desenabr@usc.edu

Copyright and License Information
=================================

Copyright (C) 2022
Andrew D. Smith and Guilherme de Sena Brandine

Authors of DNMTools: Andrew D. Smith and Guilherme de Sena Brandine

Contributors: Ben Decato, Meng Zhou, Liz Ji, Terence Li

This is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
