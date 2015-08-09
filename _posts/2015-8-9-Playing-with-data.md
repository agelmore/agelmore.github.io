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

Plot with genes per cluster as x axis and number of reads as y axis. I would guess that a low number of reads map to the smaller clusters. 

Plot x axis is gene length and y axis is read copy number. What did they do in the JC paper to account for gene length?


##To do
