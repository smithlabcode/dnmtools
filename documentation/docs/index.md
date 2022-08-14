# DNMTools

A set of tools for analyzing DNA methylation data from high-throughput
sequencing experiments, especially whole genome bisulfite sequencing
(WGBS), but also reduced representation bisulfite sequencing
(RRBS). These tools focus on overcoming the computing challenges
imposed by the scale of genome-wide DNA methylation data, which is
usually the early parts of data analysis.

Dnmtools supersedes our MethPipe software, but includes all the same
functionality (or soon will). The design of dnmtools differs, but most
commands within dnmtools are in one-to-one correspondence with some
tool in MethPipe.

At present, as we make the first release of dnmtools, the
documentation is structured very differently from that of MethPipe. We
are rebuilding that documentation here, but it might take a few
iterations to stabilize. Most of the old documentation still applies
for tasks that do not directly involve our tools (e.g., how to
manipulate the files with basic terminal commands). If you are already
a user of our tools, please check the documentation here in case
individual commands within dnmtools have changes to command line
arguments or file formats.

# Assumptions

Our tools were designed to run in a cluster computing context, with
many processing nodes available, and a job submission system like
PBS/SGE or SLURM. But now it is also possible to analyze methylomes
from most genomes (including human and mouse) in a local machine with
16 GB of RAM. Given that abismal is so fast and memory efficient, even
less may be required. Typically the data we deal with amounts to a
minimum of 100GB for a mammalian methylome at 10x coverage. This data
expands during processing before the final files are ready.
Intermediate files may cause this amount to more than double during
execution of any of our pipelines.

It is critical that users are familiar with bisulfite sequencing
experiments, especially the bisulfite conversion reaction, and how
this affects what we observe in the sequenced reads.

Users are assumed to be somewhat familar with a command line.
Installing dnmtools might require working with some environment
variables. The data that dnmtools is designed for is really big, so
understanding your computer is important. We also assume users know
what standard output and standard error are.

# Citation information

If you find our tools helpful, please cite us. Many of the commands
within dnmtools were introduced in their own publications. However,
several improvements and additions have been made over the years.
Below are the relevant citations.

**dnmtools**

Qiang Song, Benjamin Decato, Elizabeth E Hong, Meng Zhou, Fang Fang,
Jianghan Qu, Tyler Garvin, Michael Kessler, Jun Zhou, and Andrew D
Smith. *A reference methylome database and analysis pipeline to
facilitate integrative and comparative epigenomics.* PloS ONE,
8(12):e81148, 2013

**abismal**

Guilherme de Sena Brandine, and Andrew D Smith *Fast and
memory-efficient mapping of short bisulfite sequencing reads using a
two-letter alphabet.* NAR Genomics and Bioinformatics, 3(4), lqab115,
2021.

**amrfinder**

Fang Fang, Emily Hodges, Antoine Molaro, Matthew Dean, Gregory J
Hannon, and Andrew D Smith. *Genomic landscape of human
allele-specific dna methylation.* Proceedings of the National Academy
of Sciences, 109(19):7332–7337, 2012

**mlml**

Jianghan Qu, Meng Zhou, Qiang Song, Elizabeth E Hong, and Andrew D
Smith. *Mlml: Consistent simultaneous estimates of dna methylation and
hydroxymethylation.* Bioinformatics, 29(20):2645–2646, 2013

**radmeth**

Egor Dolzhenko and Andrew D Smith.  *Using beta-binomiaregression for
high-precision differen tial methylation analysis in multifactor
whole-genome bisulfite sequencing experiments.* BMC bioinformatics,
15(1):1–8, 2014

**pmd**

Benjamin E Decato, Jianghan Qu, Xiaojing Ji, Elvin Wagenblast, Simon
RV Knott, Gregory J Hannon, and Andrew D Smith. *Characterization of
universal features of partially methylated domains across tissues and
species.* Epigenetics & Chromatin, 2020

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

Essential contributors: Ben Decato, Meng Zhou, Liz Ji, Terence Li,
Jenny Qu, Qiang Song and Fang Fang

This is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
