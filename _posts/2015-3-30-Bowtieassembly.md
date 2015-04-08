---
layout: post
title:  "Bowtie assembly"
date:   2015-3-30
comments: true
---

I wasn't getting good assembly of Fusobacterium with the *de novo* assembly, so I've decided to try an assembly using a reference genome. Briefly, I will use bowtie to extract reads (from the HMP dataset) that align to Fusobacterium (from a reference database). 

I want to start by re-making my Fuso database. I'm pretty sure some draft genomes have been added since I first made the database last summer. I also want to have a pipeline so I can do this quickly in the future. All genomes are stored on the [ncbi website](ftp://ftp.ncbi.nih.gov/genomes/Bacteria).

First make database with full length Fuso genomes:

~~~~
cd 


~~~~

###Database with draft genomes added

This is a little tricker since they're a lot of them. I made a text file that has the location of all the draft genomes on the NCBI site called FASTA_location. 

{ % gist b8a98a9745241f8ab7e1 % }
~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/extract/draft

wget -i FASTA_location
~~~~