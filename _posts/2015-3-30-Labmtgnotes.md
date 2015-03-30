---
layout: post
title:  "Notes on secondary assembly"
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

Somehow when I did this the first time I was off on the cluster numbers. I must have made a bug in making the histogram, but I didn't save that code because I did it so quickly. Lesson learned, document ALL code. But this could be one explanation for why my secondary assembly didn't work! 