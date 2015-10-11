---
layout: post
title:  "Playing with data from alignment to pangenome"
date:   2015-8-5
comments: true
---

Now that I have a pipeline to align HMP reads to my pangenome, I have to normalize the data before I can figure out which genes are mapping. There are two steps to normalization: 1) normalize by gene length, 2) normalize depth to single copy core genes.

# 1) Normalize by gene length

The genes in the pangenome clusters are of varying lengths. Because of this, longer genes will map more reads (more places to map), and will therefore appear to be more common than shorter genes. Each cluster has multiple genes from the different input genomes that may be slightly different lengths (especially since some may be gene fragments). Due to this variation, I will normalize read counts to each gene individually, and then add up normalized read counts per cluster. I want to add this normalization step to my python script "sharedfile2.py" that takes my cluster index and read index to create a shared file.

The reference sequence length can be found in the same file header, or by using samtools faidx command to generate an index file. 

~~~~
#all.t0.fna is the fasta file with ALL reference genes that was used to align reads
samtools faidx all.t0.fna  

#cut out the two columns that contain reference sequence name (same as in the sam file) and sequence length
cut -f1,2 all.t0.fna.fai > all.t0.fna.cut12.fai

#Using R to merge files

x<- read.delim(file="all.t0.fna.cut12.fai", header=F)
y<- read.delim(file="all.t0.SRS013502.mapped.index", header=F)
colnames(x)<- c("seq","LN")
colnames(y)<- c("read","seq")
z<- merge(x,y,by="seq")
write.table(z, file="all.t0.SRS013502.mapped.lengths.index", quote=F, row.names=F, col.names=F, sep="\t")

~~~~

Easy. Now, using this new index file in my sharedfile.py script, I will add (1/seq length) to the count total instead of just a raw count. This will normalize each sequence for its individual sequence length. The total will be normalized for each sequence. I need to make a way to check that my script works.





