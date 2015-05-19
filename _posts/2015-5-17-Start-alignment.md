---
layout: post
title:  "Start aligning to pangenome"
date:   2015-5-16
comments: true
---

I need a pipeline to align reads to the entire set of genes and then assign each read to one of the 6500 clusters.

Outline:

1. combine genes that clustered into single file

2. create index file with gene names and cluster number

3. bwa reads to combined gene file

4. assign reads to cluster using index


##Combine genes into single file
~~~~
mkdir /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa

cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/ffn_homologues/NZACET_f0_1taxa_algOMCL_e0_/nucleotide

cat *.fna >> /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/all.t1.fna
~~~~

##Create index file
~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/ffn_homologues/NZACET_f0_1taxa_algOMCL_e0_/nucleotide

ls | sort -n > /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/cluster.names

I wonder if I can use the mothur command make.groups somehow. A group file is just a list of all the sequences and their group name. I could also easily write a script for this, probably a one-liner.
~~~~

##Run BWA

~~~~
bwa index all.t1.fna

#start with one tongue sample
wget http://downloads.hmpdacc.org/data/Illumina/tongue_dorsum/SRS013502.tar.bz2
tar jxvf SRS013502
cd SRS013502

#run bwa. without the -a option this should output only best matches.
bwa mem -M -t 16 all.t1.fna SRS013502/SRS013502.denovo_duplicates_marked.trimmed.1.fastq SRS013502/SRS013502.denovo_duplicates_marked.trimmed.2.fastq > all.t1.SRS013502.sam





