---
layout: post
title:  "Align reads to pangenome - updated"
date:   2015-7-8
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

cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/ffn_t0_homologues/NZACET_f0_0taxa_algOMCL_e0_/nucleotide

cat *.fna >> /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/all.t0.fna
~~~~

##Create index file
~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/ffn_homologues/NZACET_f0_1taxa_algOMCL_e0_/nucleotide
/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/ffn_t0_homologues/NZACET_f0_0taxa_algOMCL_e0_/nucleotide


#create a "group" file. Every line is a unique sequence and the cluster that it belongs to. I did the same for the t0 clusters
for f in *.fna; do cat $f | awk 'NR % 2 ==1' | awk -v var="$f" '{gsub(/>/,"")} {OFS="\t"} {print $1, var}' >> t0.index; done

#looks like this. The gi number specific to the sequence, the ref name of which genome (or draft) it came from, and the cluster number separated by a tab. 
head t1.index
gi|254303548|ref|ZP_04970906.1|	10000_ZP_04970906.1.fna
gi|422338993|ref|ZP_16419953.1|	10000_ZP_04970906.1.fna
gi|421527171|ref|ZP_15973775.1|	10000_ZP_04970906.1.fna
gi|254303549|ref|ZP_04970907.1|	10001_ZP_04970907.1.fna
gi|421527172|ref|ZP_15973776.1|	10001_ZP_04970907.1.fna
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
bwa mem -M -t 16 all.t0.fna ../SRS013502/SRS013502.fastq > all.t0.2.SRS013502.sam
~~~~

##Run Bowtie2

For some reason the sam file was acting weird with the bwa data so I'm trying with bowtie2 instead since I've used it before.

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0

bowtie2-build all.t0.fna all.t0.fna

bowtie2 all.t0.fna -q ../SRS013502/SRS013502.fastq -p 16 -S all.t0.SRS013502.sam  
~~~~

#Analyze alignment output

Start with the t0 output to see if reads are mapping to the singletons.

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0

#sam to bam and pull out mapped genes
samtools view -bT all.t0.fna all.t0.SRS013502.sam > all.t0.SRS013502.bam

#use the -F4 option. The -F option removes the specified FLAG. The 4 flag is unmapped reads. 
samtools view -F4 all.t0.SRS013502.bam > all.t0.SRS013502.mapped.sam

#cut out the read name and reference name
cut -f1,3 all.t0.SRS013502.mapped.sam > all.t0.SRS013502.mapped.index
~~~~

By removing unmapped reads, `427357` out of the original `33908624` remained. This means about `1.3%` of the reads in the sample mapped to fuso sequences. I was getting more than that before, but this is only one sample where I was using a pool of samples before. We expect there to be variation among samples. 


This index file contains the reads mapping to reference sequences. Together with the t0.index file containing the reference sequences in each cluster, I can create a script that will count the reads per cluster, or a "shared file".


##Creating shared file

I wrote a python script which counts the number of reads that map to each cluster in the pangenome. The input files are an index of reads mapping to sequences and an index of which sequences are in each cluster. 

The script is in a repository on my personal github [here](https://github.com/agelmore/Pangenome/blob/master/sharedfile.py).

The repository is saved on axiom here: /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/Pangenome

~~~~
#run sharedfile.py. takes about 5 minutes with a single sample

python2.7 sharedfile.py /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/t0.index /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/all.t0.SRS013502.mapped.index  /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/all.t0.SRS013502.mapped.out
~~~~


