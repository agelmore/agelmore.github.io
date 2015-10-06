---
layout: post
title:  "Playing with data from alignment to pangenome"
date:   2015-8-5
comments: true
---

In my previous posts I did pangenome construction, HMP read alignment to pangenome, and organized output into a "shared file" (with the script [sharedfile.py](https://github.com/agelmore/Pangenome/blob/master/sharedfile.py)). Here, I want to play with this output to understand my data better and come up with a plan for comparing samples.

##Generate shared file

Here's a summary of what I have done up to this point. See [this post](http://agelmore.github.io/2015/07/08/Start-alignment.html) for the alignment and data organization pipeline. 

#Files

These are the most important files in the analysis.

File | suffix | description | column 1 | column 2 | column 3
:---------------|:--------:|:--------:|:--------:|:--------:|:--------:
Cluster index | *.index | genes matched with cluster from pangenome | gene name | cluster name 
Read index | *.mapped.index | reads matched with gene from sam file | read name | gene name 
Shared file | *.mapped.out | Count of reads and genes in each cluster | cluster name | read count | gene count 




##Read copy number per cluster

Now I can start to look at the distribution of reads on different gene clusters using the shared file I created with sharedfile.py. 

~~~~
setwd("~/Documents/Schloss/Fuso/Pangenome/")

x<- read.delim(file="all.t0.SRS013502.mapped.out", header=T)
plot(x$genecount[x$genecount<6], x$readcount[x$genecount<6], type="p", main="Do reads map to clusters with few genes?", xlab="number of gene copies", ylab='number of reads mapping to cluster')
~~~~ 

![Reads mapping to rare clusters]({{ site.url }}/images/genesvsreads5.png)

Plot x axis is gene length and y axis is read copy number. What did they do in the JC paper to account for gene length?

Kathryn multiplies by 100 (read length) and divides by the gene length. I could somehow average the gene length.

Figure out if they're real genes in the single gene clusters.

##Pangenome rarefaction

We are curious to know more about the pangenome. There are so many rare genes, so we want to know if those are evenly spread through the input genomes. We can learn about the distribution of rare genes by using rarefaction. The technique is traditionally used in ecology to calculate species richness. It creates a curve with the number of species (in this case gene clusters) as a function of the number of samples (in this case genomes in pangenome). A sharp slope indicates high diversity because adding few genomes drastically increases the number of clusters; a flat line would happen when adding new genomes to the pangenome does not contribute any rare genes (all of the novel genes have been discovered).

To create this curve, I will use the mothur rarefaction.shared command which only requires a shared file. In this file I will have gene clusters in place of OTUs (columns) and genome names in place of sample names (rows). I will have to reorganize the data in the file `/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/ffn_t0_homologues/NZACET_f0_0taxa_algOMCL_e0_.cluster_list` to make this pangenome shared file. I can do this using a script from the get_homologues package.

~~~~
compare_clusters.pl -o faa_compare2 -d ffn_t0_homologues/NZACET_f0_0taxa_algOMCL_e0_ -t 0 -m
~~~~

The output file `/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_compare/pangenome_matrix_t0.tab` is basically this shared file, but without the "label" and "numOtus" column. These are required for the mothur rarefaction.shared command, so I added these manually.

Run mothur:

~~~~
rarefaction.shared(shared=pangenome_matrix_t0_manedits.tab, label=0.00)
~~~~

Plot in R:

~~~~
setwd("~/Documents/Schloss/Fuso/Pangenome/analysis")
x<- read.delim(file="pangenome_matrix_t0_manedits.shared.rarefaction", header=T)
plot(x$numsampled,x$X, xlab="Number of genomes", ylab="Number of clusters", main="Rarefaction curve with t=0 clustering, singletons included")
~~~~

![Rarefaction curve - t=0, singletons]({{ site.url }}/images/rarefaction.t0.singletons.png)

The line doesn't flatten out which means that each genome is adding a new set of novel genes to the pangenome. This shouldn't be the case. I'm going to try making curves of the pangenome clustered at t=1.

###Rarefaction t=1

~~~~
compare_clusters.pl -o faa_compare_t1 -d faa_draft_t1/NZACET_f0_1taxa_algOMCL_e0_ -t 0 -m 
~~~~

Manually edit to add label column and changed file name to `pangenome_matrix_t1_manedits.tab`

Mothur:

~~~~
rarefaction.shared(shared=pangenome_matrix_t1_manedits.tab, label=0.00)
~~~~

Add to R plot:

![Rarefaction curve - t=1, no singletons]({{ site.url }}/images/rarefaction.t1.png)
