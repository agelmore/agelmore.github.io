---
layout: post
title:  "Automating alignment pipeline"
date:   2015-10-18
comments: true
---

I now have a pipeline to create a gene abundance profile to the Fuso pangenome in one sample. I would like to make this more streamlined so I can do it in many samples in various body sites and compare.

#Summary of pipeline with sample2

##Files

What files do I need that will be the same for all samples?

All genes file: `/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/all.t0.fna`
#This is a fasta file that has all the genes from the 33 fuso genomes combined into a single file
Cluster index file: `/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/t0.index'
#Index file with sequence name in first column and cluster name in second column

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

##Bowtie2

~~~~
cd ..
mkdir bowtie; cd bowtie

bowtie2-build ../panfiles/all.t0.fna all.t0.fna

bowtie2 all.t0.fna -q ../data/SRS013705/SRS013705.fastq -p 16 -S all.t0.SRS013705.sam 
~~~~

##Create shared file

**Put summary of shared file here**

~~~~
cd ..
mkdir shared; cd shared

#sam to bam and pull out mapped genes
samtools view -bT ../panfiles/all.t0.fna ../bowtie/all.t0.SRS013705.sam > all.t0.SRS013705.bam

#use the -F4 option. The -F option removes the specified FLAG. The 4 flag is unmapped reads. 
samtools view -F4 all.t0.SRS013705.bam > all.t0.SRS013705.mapped.sam

#remove bam file to save space
rm all.t0.SRS013705.bam

#cut out the read name and reference name
cut -f1,3 all.t0.SRS013705.mapped.sam > all.t0.SRS013705.mapped.index

#use the faidx command to get sequence lengths from the reference file
samtools faidx ../panfiles/all.t0.fna  

#cut out the two columns that contain reference sequence name (same as in the sam file) and sequence length
cut -f1,2 ../panfiles/all.t0.fna.fai > all.t0.fna.lengths

#Using R to merge files

x<- read.delim(file="all.t0.fna.lengths", header=F)
y<- read.delim(file="all.t0.SRS013705.mapped.index", header=F)
colnames(x)<- c("seq","LN")
colnames(y)<- c("read","seq")
z<- merge(x,y,by="seq")
write.table(z, file="all.t0.SRS013705.mapped.lengths.index", quote=F, row.names=F, col.names=F, sep="\t")

~~~~

Run shared file script from sample2 folder:

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
hist(x[x$readcount>1,'readcount'], xlim=c(1,40), breaks=5000,xlab="Gene coverage (per base)", main="Frequency of gene abundance, normalized by gene length, SRS013705", file ="abundance.SRS013705.png')
dev.off()

~~~~


![Histogram of gene abundance SRS013705]({{ site.url }}/images/abundance.SRS013705.png)











