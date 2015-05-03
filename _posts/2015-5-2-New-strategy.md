---
layout: post
title:  "New strategy"
date:   2015-5-2
comments: true
---

I gave chalk talk for lab meeting last Tuesday, and we decided I should take a new direction with my project, outlined here. 

#Assembly-based strategy

Until this point I've been pursuing an assembly strategy. I hoped that I would be able to assemble the metagenomic reads well enough that I could separate contigs into different strains. Once I had created fully assembled genomes, I would be able to compare genes and abundance at different body sites by mapping to the strains that I had assembled.  I spent a lot of time learning about different assemblers, normalization, de novo assembly, and reference-based assembly. In the end, I created my best assembly using a reference based approach: bowtie extracted reads, digital normalization, megahit assembly. 

The problem is that even my best effort at an assembly really isn't that good. Problems with the assembly based approach:

1. **There aren't enough reads in each sample to get a full length assembly.** Despite pooling 20 samples, I still have huge gaps. I need more information if I'm going to be able to generate long enough contigs to separate out strains using an assembler.

2. **By assembling reads into contigs, I lose variation in the abundance of different genes.** A gene variant with 100x read coverage and a gene variant with 10x read coverage will both collapse down to a single contig. In order to understand the variation, I will have to go back and look at the reads anyways. Therefore, assembling is an unnecessary step.

3. **Some strains may have regions of the genome that are identical.** When assembling, reads from different strains may be combined into a single assembled contig(s). When I use CONCOCT or other software (ESOM) that can separate the assembled contigs, contigs that are shared between strains will be clustered together. Therefore, sequence similarity will make it impossible to use this method to separate strains and assemble complete genomes. 

#Pangenome read variant calling

My new strategy is to look at variation in Fusobacterium genes using abundance of read variants. 

Features of this strategy:

1. **Creation of a Fusobacterium Pangenome**. NCBI has 4 full Fusobacterium genomes and many draft genomes. In this dataset there are thousands of Fuso genes that span many species and strains. I will create a single database that is a unique list of the Fuso genes present in any strain.

2. **Map reads to Pangenome**. Using Bowtie2.

3. **Read variant calling**. Using samtools variant calling, I will see which genes have the highest read variation. I will also look at read coverage on each gene. Low variation genes with high abundance are genes that are identical between strains/species. Low variation genes with low abundance may be unique genes found in only some strains/species. We are most interested in high variation genes with high abundance because these are genes that are probably under selection pressure in different strains. 

4. **Gene variation across body sites, individuals, and disease states**. Once I have developed a pipeline for looking at variation in the pangenome in different samples, I will be able to quickly run different samples to compare gene variation.


###Possible issues

As I'm brainstorming I thought of a few possible issues with my new strategy.

1. **How do I make sure the mapped reads aren't other bugs (not Fuso) with similar sequences?** If my Fuso pangenome contains genes that are very similar to genes in other species, I might find those genes to have artificially high read abundance/variation. In a recent paper by [Greenblum et al.](http://www-ncbi-nlm-nih-gov.proxy.lib.umich.edu/pubmed/25640238), reads were aligned to reference genomes to survey gene copy number. They performed simulation-based analyses to validate that reads originating from a certain gene mapped to the correct genome and gene. The simulated 75bp reads and mapped them to their 260 reference genome database (which had been clustered into genome clusters based on sequence similarity of marker genes). They got an error rate of about 11.8% (including synthetic reads that mapped to the wrong gene/genome or didn't map at all). They greatly decreased this error rate by removing some reference genomes from the database. They also identified the optimal maximum edit distance (MED) in the alignment at 5. A more strict alignment excluded reads with sequence errors or small amounts of variation; a more lenient alignment increased false positives. **I will need to do some sort of simulation analysis to find the best MED and identify genes that have a high false positive rate to possible remove from the analysis.**


2. **Novel genes not present in database**. There is a possibility that there are Fusobacterium reads in my samples that correspond to novel genes that have not been sequenced. Unfortunately, with this pipeline those genes will be overlooked. Assembly would be the only way to find novel genes. Because it is much easier to assemble single genome sequences, in the future I can sequence cultured isolates in the lab to find new genes. I can add these to my pangenome. 









