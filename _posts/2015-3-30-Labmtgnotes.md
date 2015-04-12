---
layout: post
title:  "Redo secondary assembly"
date:   2015-3-29
comments: true
---

I got a lot of good comments during my code review in lab meeting. One was regarding the reason I picked cluster 17 and 21 for the secondary assembly. 

In a [previous post](http://agelmore.github.io/2015/02/16/Concoct-1kb.html) where I initially ran through the CONCOCT pipeline on this dataset, I created a histogram showing the frequency of Fuso in each of the clusters. This graph:

![Cluster Histogram]({{ site.url }}/images/clusterhistogram.png)

But Niel pointed out that frequency isn't the best way to show this graph, because there are probably different numbers of total contigs in each bin. I recreated this graph in R with the percent of Fuso contigs instead of frequency.

~~~~
library(plyr)
y<- read.table("blast/blast.concoct.merged.1kb.txt", header=T)  #contigs with cluster that blasted to fuso
x<- read.delim("concoct/1kb/concoct-output/clustering_gt1000.csv", sep=",", header=F) #contigs with cluster total

fusocounts<-count(y$cluster) #use the plyr package
colnames(fusocounts)<- c("cluster","fuso")
totalcounts<- count(x$V2)
colnames(totalcounts)<- c("cluster", "total")
mer<- merge(totalcounts,fusocounts,by="cluster")    #table with total contigs and fuso contigs per cluster
mer$percent <- mer$fuso/mer$total *100
png('plot2.png')
plot(mer$cluster,mer$percent,type="h", lwd=10, main="Percent of Blasted Fuso contigs in CONCOCT clusters", xlab="Bin number", ylab="Percent contigs containing Fuso", las=2)
png(percentplot, file="plot.png")
dev.off()
~~~~ 

Here is the graph:

![Cluster Histogram]({{ site.url }}/images/plot2.png)

This graph shows that even when normalized for the number of contigs, clusters **18 and 22** have the most Fuso. 

Somehow when I did this the first time I was off on the cluster numbers. I must have made a bug in making the histogram, but I didn't save that code because I did it so quickly. Lesson learned, document ALL code. But this could be one explanation for why my secondary assembly didn't work! Since I have all **that** code documented, it's easy to do again.


~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/run2/concoct/1kb/assembly2
awk -F , '$2 == "18"' ../concoct-output/clustering_gt1000.csv > cluster.18 
cut -d "," -f1 cluster.18 > cluster.18.cut
awk -F , '$2 == "22"' ../concoct-output/clustering_gt1000.csv > cluster.22
cut -d "," -f1 cluster.22 > cluster.22.cut
cat cluster.18.cut cluster.22.cut > cluster.18.22

mothur '#get.seqs(accnos=cluster.18.22, fasta=megahit.1000.contigs_c10K.fa)'
mv megahit.1000.contigs_c10K.pick.fa cluster.18.22.fa

bowtie2-build cluster.18.22.fa cluster.18.22
bowtie2 cluster.18.22 -q $HMP/D1.tongue/run2/cat/All.D1.tongue.run2.cat.fq -p 16 -S cluster.18.22.sam 

#make BAM file
samtools view -bT cluster.18.22 cluster.18.22.sam > cluster.18.22.bam

#use the -F4 option. The -F option removes the specified FLAG. The 4 flag is unmapped reads. 
samtools view -F4 cluster.18.22.bam > cluster.18.22.mapped.sam

#cut out name, sequence, and quality from sam file
cut -f1,10,11 cluster.18.22.mapped.sam > cluster.18.22.mapped.cut.sam

awk '{print "@"$1"\n"$2"\n""\+""\n"$3}' cluster.18.22.mapped.cut.sam > cluster.18.22.mapped.cut.fastq

khmerEnv
cd $HMP/D1.tongue/run2/concoct/1kb/assembly2
mkdir DN.18.22
cd DN.18.22
python2.7 /mnt/EXT/Schloss-data/amanda/Fuso/khmer/khmerEnv/bin/normalize-by-median.py -C 20 -k 21 -x 1e9 ../cluster.18.22.mapped.cut.fastq -s cluster1822.D1.Tongue.run2.savetable -o cluster.18.22.mapped.cut.normalized.fastq

cd /mnt/EXT/Schloss-data/amanda/Fuso/megahit/megahit

python ./megahit -m 45e9 -r $HMP/D1.tongue/run2/concoct/1kb/assembly2/DN.18.22/cluster.18.22.mapped.cut.normalized.fastq --cpu-only -l 101 -o $HMP/D1.tongue/run2/concoct/1kb/assembly2/megahit/normalized/18.22.2

~~~~


It keeps failing on the assembly. No idea why. Here is the tail of the outfile:

~~~~
[Sun Apr  5 15:08:23 2015]: Extracting iterative edges from k = 41 to 51
[Sun Apr  5 15:08:23 2015]: Building graph for k = 51
[Sun Apr  5 15:08:23 2015]: Assembling contigs from SdBG for k = 51
[Sun Apr  5 15:08:24 2015]: Extracting iterative edges from k = 51 to 61
[Sun Apr  5 15:08:24 2015]: Building graph for k = 61
[Sun Apr  5 15:08:24 2015]: Assembling contigs from SdBG for k = 61
[Sun Apr  5 15:08:24 2015]: Extracting iterative edges from k = 61 to 71
[Sun Apr  5 15:08:24 2015]: Building graph for k = 71
[Sun Apr  5 15:08:24 2015]: Assembling contigs from SdBG for k = 71
[Sun Apr  5 15:08:24 2015]: Extracting iterative edges from k = 71 to 79
[Sun Apr  5 15:08:24 2015]: Building graph for k = 79
Error occurs when running "builder build" for k = 79
[Exit code -8]
Number of CPU threads 16
qsub working directory absolute is
/mnt/EXT/Schloss-data/amanda/Fuso/megahit/megahit

~~~~