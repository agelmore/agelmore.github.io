---
layout: post
title:  "Blasting for Fuso"
date:   2015-1-28
comments: true
---

While I'm waiting for the CONOCOCT pipeline to run on axiom, I thought I would do a couple BLAST alignments of the megahit assembly to some Fuso genomes that I downloaded a while ago.

The NCBI has over 1000 full bacteria genomes that have been properly annotated. They can be downloaded FTP from ftp://ftp.ncbi.nih.gov/genomes/Bacteria/. I pulled out all of the 4 Fuso genomes as well as all of the Fuso draft sequences. 

First I combined all of these references into a single file. Then I made the reference genome into a database and then blast with contigs as query. I used the option "max_target_seqs" so that each contig would only be mapped once. This will avoid contigs that map to multiple Fuso genomes (such as 16S):

~~~~
makeblastdb -in Fuso.all.db -dbtype nucl -out Fuso.all.db.make

blastn -db Fuso.all.db.make -query /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/run2/megahit/megahit.contigs_c10K.fa -out blast.fuso.all -evalue 1e-5 -outfmt 6 -num_threads 16 -max_target_seqs 1
~~~~

There were **40015** contigs that matched above the 1e-5 evalue cutoff. That is **2.85%** of the 1403622 contigs. 
