---
layout: post
title:  "Normalizing read counts"
date:   2015-10-11
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

Easy. Now, using this new index file in my sharedfile.py script, I will divide the read count of each sequence by its length and multiply by 100. This will normalize each sequence so that the totals are a per base coverage. I actually decided to multiple the total decimal by 100 to save computational time (you can pull out the 100 in addition). I made the changes to the [sharedfile2.py](https://github.com/agelmore/Pangenome/blob/master/sharedfile2.py) script that are saved on my github page. 

~~~~
python2.7 sharedfile2.py /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/t0.index /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/all.t0.SRS013502.mapped.lengths.index /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/all.t0.SRS013502.mapped.out3

~~~~

It takes a little longer than before because it has to compute the float division in each loop. Oh well, I just won't tell Kathy. It works!

And a little plot to look at gene coverage of the normalized data:

~~~~
setwd("~/Documents/Schloss/Fuso/Pangenome/")

x<- read.delim(file="all.t0.SRS013502.mapped.out3", header=T)

hist(x[x$readcount>1,'readcount'], xlim=c(1,40), breaks=5000,xlab="Gene coverage (per base)", main="Frequency of gene abundance, normalized by gene length")

~~~~

![Histogram of gene abundance]({{ site.url }}/images/freq.normalized.png)

As we hoped to see, the frequency of the different gene clusters are pretty normally distributed. I normalized to gene length and multiplied by read lenth, so the y axis is the single base coverage of Fuso genes in this sample.  I have not yet normalized to single copy core genes, so the center of the distribution (around 17) doesn't mean much. 

*The spike in clusters that are in low abundance may just be sequencing errors or misalignments?*

The clusters that are in high abundance may indicate genes that are present in all Fuso strains, and may be under selective pressures. It would be interesting to look at what those genes are. Next, I will do the same pipeline for other samples to see if the same genes are present at high abundance in different individuals and body sites. 

# 2) Normalize by single copy core genes

The CONCOCT software which I used last fall has a pipeline to validate purity of clusters using single copy core genes. The program uses a list of 36 COGs that are present only one time in at least 97% of the genomes on NCBI. I figured that these would be good candidates for normalization as the copy number of these genes should indicate the number of strains in the sample. 

I'm going to try to adapt the scripts from the CONCOCT pipeline to create a table with the number of COGs in my pangenome. 

~~~~
#have to change the fasta names to only include gi reference name and not the comments
awk '{if(NR%2==1){print $1"_"}else{print $1}}' all.t0.faa > all.t0.rename.faa

#alter names so the CONCOCT script works
sed -i 's/|/_/g' all.t0.rename.faa 
sed -i 's/|/_/g' t0.index.csv

$CONCOCT/scripts/RPSBLAST.sh -f all.t0.rename.faa -p -c 8 -r 1

python $CONCOCT/scripts/COG_table.py -b all.out -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/bwa/t0/t0.index.csv --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > cogtable_scg.tsv

Rscript $CONCOCT/scripts/COGPlot.R -s cogtable_scg.tsv -o cogtable_scg.pdf
~~~~


![Pangenome clusters containing single copy core genes]({{ site.url }}/images/cogtable_scg.png)


Awesome!!! It looks like there are 35 clusters that have 30 or more copies of different single copy core genes. There is only one gene (COG0504 - CTP synthase (UTP-ammonia lyase)) that is not represented in >30 copies. **I can use this information for normalization.** The clusters that contain single copy core genes should be present in the sample one time for every strain. If I divide the abundance of all genes by the abundance of these clusters, I will get an idea of the frequency of rare gene abundance in different samples. *Hopefully, these clusters are present at an even abundance in the sample*. 

This is also a good validation of my input genomes. I had 33 genomes that I used to create the pangenome, so clusters with ~33 copies of the single copy genes are represented in the majority of my input genomes. This suggests that even the draft genomes are complete without many repeat genes. 

*What does this say about the pangenome clusters?* There are 14 clusters with a few copies of the COGs. These are repeats from clusters which have 30-33 copies. The COGs should only be represented on time in every genome, so if they are hitting more than one cluster, there is something wrong with my pangenome. Could this be a misalignment in the CONCOCT pipeline? Errors in the input genomes?


###Normalize sample using SCG clusters

Extract cluster names from COG table using R. Adapted from COGplot.R script.

~~~~
tab<- read.delim(file="cogtable_scg.tsv",row.names=1, header=TRUE)
ecogs <- tab[,3:ncol(tab)]
maxecogs <- max(ecogs)
maxecogs <- max(ecogs)
sumCogs <- rowSums(ecogs)
ecogs$sum <- sumCogs
ecogs <- subset(ecogs,ecogs$sum > 25) #only want high copy number clusters
ecogs.order <- ecogs[order(ecogs$sum),]
names <- row.names(ecogs.order)
write.table(names, file="clusternames", row.names=FALSE, quote=FALSE)
~~~~

Merge SCG clusters with abundance table to get abundance of SCG clusters.

~~~~
setwd("~/Documents/Schloss/Fuso/Pangenome/analysis")
x<- read.delim(file="pangenome_matrix_t0_manedits.shared.rarefaction", header=T)
z<- read.delim(file="clusternames", header=TRUE)
colnames(z)<- ("clustername")
merg<- merge(x,z,by="clustername")
mean(merg$readcount)

~~~~

`21.05792` is the average abundance of the SCG clusters. I will use this number to normalize this sample. This number will be different for each sample I do, using the average abundance of these clusters.

*Why is there variation in these clusters across the sample? Sequencing error? Alignment error?*



