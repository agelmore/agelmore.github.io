---
layout: post
title:  "Creating a pangenome"
date:   2015-5-3
comments: true
---


A pangenome of a species is the sum of the core genes found in all strains and strain-specific genes found in some but not all strains [(review)](http://www-ncbi-nlm-nih-gov.proxy.lib.umich.edu/pubmed/19086349).

I will be using the software package [GET_HOMOLOGUES](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3837814/) (Bruno Contreras-Moreira and Vinuesa, 2013. *Genome Research*). 

![GET_HOMOLOGUES Flowchart ]({{ site.url }}/images/get_homologues_flowchart.png)

The software package is open-source, automatic, customizable, and can handle large jobs. All the scripts are written in perl, yay. Basically, genes input in FASTA format are BLASTed against each other. Based on BLAST hits, genes are clustered into orthologous groups using one of three popular clustering algorithms: COGtriangles, Bidirectional Best hit, and OrthoMCL (recommended by Evan Snitkin). The final result is a fasta file with all the pangenomic genes (other formats and statistics optional). The paper also showed that the pipeline works for draft genomes, which is good because I have a lot of them.

##Install

~~~~
cd /mnt/EXT/Schloss-data/amanda/get_homologues

wget http://maya.ccg.unam.mx/soft/get_homologues-x86_64-20150306.tgz
tar xvfz get_homologues-x86_64-20150306.tgz
./install.pl

#test it
./get_homologues.pl -d sample_buch_fasta 
~~~~

##Format FASTA

The best way to input files into GET_HOMOLOGUES is by putting all of the sequences in a single folder in FASTA format. Each genome is in a separate file. This way, if more genomes are added later (or if we sequence some in lab), they can be easily added to the pangenome. 

I'm going to start with a pangenome of just F. nucleatum using complete genomes. See if that works and I can add the draft genomes later. I'm downloading all the ncbi entries into a single folder.

I tried to run the program using the fasta files from ncbi, but those haven't been annotated like the faa (amino acid) files have. I'm going to try running get_homologues with the amino acid files. This will produce a protein sequences pangenome, but maybe I can translate it back later? 

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/fasta

wget -i complete_nucleatum_faa.txt
~~~~

##Run

~~~~
get_homologues.pl -i faa
~~~~








