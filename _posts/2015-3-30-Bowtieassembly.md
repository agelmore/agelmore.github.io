---
layout: post
title:  "Bowtie assembly"
date:   2015-3-30
comments: true
---

I wasn't getting good assembly of Fusobacterium with the *de novo* assembly, so I've decided to try an assembly using a reference genome. Briefly, I will use bowtie to extract reads (from the HMP dataset) that align to Fusobacterium (from a reference database). 

I want to start by re-making my Fuso database. I'm pretty sure some draft genomes have been added since I first made the database last summer. I also want to have a pipeline so I can do this quickly in the future. All genomes are stored on the [ncbi website](ftp://ftp.ncbi.nih.gov/genomes/Bacteria).

###First make database with full length Fuso genomes

I downloaded all of the full length genomes from ncbi since they don't take up too much space. I combined all the fuso genomes (and full length plasmids) into a single file. Here are the sequence names:

{% gist 7cdeb4fc470cb793259c %}

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/extract

wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/*    #download all genomes
for f in Fuso*; do cat $f/*.fna >> Database/fusodb.fa; done

~~~~

###Database with draft genomes added

This is a little tricker since they're a lot of them. I made a text file that has the location of all the draft genomes on the NCBI site called FASTA_location. Each of these genomes has multiple scaffolds because they aren't complete genomes. I combined all of them into a single file. 

{% gist b8a98a9745241f8ab7e1 %}

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/extract/draft

wget -i FASTA_location
cat *.fna >> draftdb.fa
~~~~

And then combine the databases into one...

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database
cat fusodb.fa ../draft/draftdb.fa > fusodb.complete.fa
~~~~







