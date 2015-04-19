---
layout: post
title:  "Coverage map on F. nucleatum"
date:   2015-4-19
comments: true
---

Post to make a plot of where HMP reads map on F. nucleatum genome. This will give us an idea of where the reads are aligning to the reference genomes and what coverage level we have. 

I'm going to align the HMP sample reads to a single F. nucleatum genome (already have the database built here `/mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single`). 

Run blast:

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/coverage

blastn -db /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single -query $HMP/D1.tongue/run2/cat/All.D1.tongue.run2.cat.fa -out blast.single -evalue 1e-5 -outfmt 6 -num_threads 16 -max_target_seqs 1
~~~~

Or bowtie:

~~~~
bowtie2 /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single -q /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/run2/DN/All.D1.Tongue.run2.norm.fq -p 16 -S bowtie.single.sam
~~~~

genomeCoveragebed is a bedtool that will calculate the coverage at each base of the designated genome. 

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/coverage

#create bam file
samtools view -Sb bowtie.single.sam > bowtie.single.bam

#format the reference genome
samtools faidx /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single.fna  

#sort the bam file
samtools sort bowtie.single.bam bowtie.single.sorted.bam

#run genomeCoverageBed. The bg option will produce a BEDGRAPH output file where the coverage of each position is listed. 
bedtools genomecov -bg -i bowtie.single.sorted.bam -g /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single.fna.fai > bowtie.single.bg

~~~~

It's all running on axiom. When it's finished I'll figure out how to plot the bedgraph output in R.