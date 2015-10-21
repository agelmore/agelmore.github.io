---
layout: post
title:  "Automating alignment pipeline"
date:   2015-10-18
comments: true
---

I now have a pipeline to create a gene abundance profile to the Fuso pangenome in one sample. I would like to make this more streamlined so I can do it in many samples in various body sites and compare.

#Summary of pipeline with sample2

1. **Creation of a Fusobacterium Pangenome**. NCBI has 3 full Fusobacterium genomes and 30 draft genomes. In this dataset there are thousands of Fuso genes that span many species and strains. I created a single database that is a unique list of the Fuso genes present in any strain.

2. **Map reads from the HMP to Pangenome**. Using Bowtie2, I will map samples one at a time to create a gene abundance profile for each individual at each body site.

3. **Normalize each sample**

4. **Create shared file with all samples** Like a mothur shared file. "OTUs" are pangenome clusters. Each sample will have a list of abundance of each cluster in units of depth per base pair in the gene.

~~~~
clustername	sample1	sample2
25519_ZP_10973098.1.fna	0.00000000	0.00000000
8061_ZP_00144933.1.fna	0.00000000	0.00000000
60135_ZP_16967793.1.fna	0.17482517	1.60390516
17632_ZP_07924364.1.fna	0.00000000	0.00000000
66745_ZP_09587444.1.fna	0.00000000	0.0000000
9899_ZP_04970805.1.fna	2.15818690	44.87292999
22298_ZP_08691379.1.fna	0.00000000	10.16260163
62861_ZP_09177555.1.fna	0.00000000	0.00000000
7195_ZP_00144067.1.fna	5.54089710	0.00000000
34487_ZP_16396549.1.fna	0.00000000	10.16260163
62846_ZP_09177540.1.fna	2.66666667	0.00000000
59236_ZP_16966894.1.fna	1.17647059	1.17647059
77016_ZP_15973812.1.fna	3.08123249	1.17647059

~~~~

##Files

What files do I need that will be the same for all samples?

All genes file: `/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/all.t0.fna`
#This is a fasta file that has all the genes from the 33 fuso genomes combined into a single file
Cluster index file: `/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/t0.index`
#Index file of pangenome with sequence name in first column and cluster name in second column

Create a new directory for the new sample

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/
mkdir sample2; cd sample2
mkdir panfiles; cd panfiles
cp /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/all.t0.fna .
cp /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/t0.index .
~~~~


##Download sample

~~~~
#Download sample
mkdir data; cd data
wget http://downloads.hmpdacc.org/data/Illumina/tongue_dorsum/SRS013705.tar.bz2
tar jxvf SRS013705.tar.bz2
cd SRS013705
cat *.1.fastq *.2.fastq > SRS013705.fastq
~~~~

**Can I make this into a bash script that inputs sample name?**

##Bowtie2 alignment

Align reads to complete gene database from all 33 genomes.

~~~~
cd ..
mkdir bowtie; cd bowtie

bowtie2-build ../panfiles/all.t0.fna all.t0.fna

bowtie2 all.t0.fna -q ../data/SRS013705/SRS013705.fastq -p 16 -S all.t0.SRS013705.sam 
~~~~

##Create read abundance profile


~~~~
cd ..
mkdir shared; cd shared

#sam to bam 
samtools view -bT ../panfiles/all.t0.fna ../bowtie/all.t0.SRS013705.sam > all.t0.SRS013705.bam

#use the -F4 option to pull out mapped reads. The -F option removes the specified FLAG. The 4 flag is unmapped reads. 
samtools view -F4 all.t0.SRS013705.bam > all.t0.SRS013705.mapped.sam

#remove bam file to save space
rm all.t0.SRS013705.bam

#cut out the read name and reference name to create index of mapped reads
cut -f1,3 all.t0.SRS013705.mapped.sam > all.t0.SRS013705.mapped.index

#use the faidx command to get sequence lengths from the reference file (to use for normalization)
samtools faidx ../panfiles/all.t0.fna  

#cut out the two columns that contain reference sequence name (same as in the sam file) and sequence length
cut -f1,2 ../panfiles/all.t0.fna.fai > all.t0.fna.lengths

#Using R to merge gene lengths with mapped reads

x<- read.delim(file="all.t0.fna.lengths", header=F)
y<- read.delim(file="all.t0.SRS013705.mapped.index", header=F)
colnames(x)<- c("seq","LN")
colnames(y)<- c("read","seq")
z<- merge(x,y,by="seq")
write.table(z, file="all.t0.SRS013705.mapped.lengths.index", quote=F, row.names=F, col.names=F, sep="\t")

~~~~

Run [shared file script]((https://github.com/agelmore/Pangenome/blob/master/sharedfile.py) from sample2 folder:

~~~~
cd ..
mkdir output; cd output
python2.7 /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/Pangenome/sharedfile2.py /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/sample2/panfiles/t0.index /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/sample2/shared/all.t0.SRS013705.mapped.lengths.index /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/sample2/output/all.t0.SRS013705.shared

~~~~

**Is there a way to call this script without the long directories?**

##Look at distribution

**R**

~~~~
cd output
png('abundance.SRS013705.png')
x<- read.delim(file="all.t0.SRS013705.shared", header=T)
hist(x[x$readcount>1,'readcount'], xlim=c(1,100), breaks=5000,xlab="Gene coverage (per base)", main="Frequency of gene abundance, normalized by gene length, SRS013705", file ="abundance.SRS013705.png")
dev.off()

~~~~


![Histogram of gene abundance SRS013705]({{ site.url }}/images/abundance.SRS013705.png)










