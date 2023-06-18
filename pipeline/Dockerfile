# syntax=docker/dockerfile:1

## Copy this file to where the directory where you have the dnmtools binary. Make sure it was built on linux x86_64.

FROM ubuntu:22.04

# install pipeline dependencies
RUN apt-get update
RUN apt-get install -y libgsl-dev libhts-dev libgomp1 samtools libcurl4 trim-galore sra-toolkit rsync
RUN rm -rf /var/lib/apt/lists/*
RUN rsync -a hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/bedToBigBed /usr/bin
RUN rsync -a hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/wigToBigWig /usr/bin

# install dnmtools and it must be build for Ubuntu
COPY dnmtools /usr/bin
