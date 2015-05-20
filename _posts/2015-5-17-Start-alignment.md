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
/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/ffn_t0_homologues/NZACET_f0_0taxa_algOMCL_e0_/nucleotide


#create a "group" file. Every line is a unique sequence and the cluster that it belongs to. I did the same for the t0 clusters
for f in *1.fna; do cat $f | awk 'NR % 2 ==1' | awk -v var="$f" '{gsub(/>/,"")} {OFS="\t"} {print $1, var}' >> t1.index; done

~~~~

##Run BWA

~~~~
bwa index all.t1.fna

#start with one tongue sample
wget http://downloads.hmpdacc.org/data/Illumina/tongue_dorsum/SRS013502.tar.bz2
tar jxvf SRS013502
cd SRS013502
cat *.1.fastq *.2.fastq > SRS013502.fastq

#run bwa. without the -a option this should output only best matches.
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t1
bwa mem -M -t 16 all.t1.fna ../SRS013502/SRS013502.fastq > all.t1.SRS013502.sam

#run again against database including singletons
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0
bwa mem -M -t 16 all.t0.fna ../SRS013502/SRS013502.fastq > all.t0.SRS013502.sam
~~~~

#Analyze BWA output

Start with the t0 output to see if reads are mapping to the singletons.

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0

#sam to bam
samtools view -Sb all.t0.SRS013502.sam > all.t0.SRS013502.bam

~~~~


ideas to use:
coverageBed -abam all.t1.SRS013502.bam -b all.t1.fna
concoct scripts - concoct basically does what I want to do. They map all the reads to the contigs and then generate a table with the coverage per contig. 
mothur - takes an group file and sequences and generates a shared file


