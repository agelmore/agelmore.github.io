---
layout: post
title:  "Project outline"
date:   2015-1-11
comments: true
---

Brief background and basic outline for my project studying *Fusobacterium nucleatum*.

##Background

In our lab, we are interested in the role of the gut microbiome in the development and progression of colorectal cancer (CRC). Patients with CRC appear to have structurally distinct gut microbiomes; for example, these patients frequently harbor an enrichment of the commensal *Fusobacterium nucleatum* along with other shifts in structure ([Kostic et al., 2013](http://www.ncbi.nlm.nih.gov/pubmed/23954159), [Kostic et al., 2012](http://www.ncbi.nlm.nih.gov/pubmed/?term=Genomic+analysis+identifies+association+of+Fusobacterium+with+colorectal+carcinoma), [Castellarin et al., 2012](http://www.ncbi.nlm.nih.gov/pubmed/?term=Fusobacterium+nucleatum+infection+is+prevalent+in+human+colorectal+carcinoma)). *F. nucleatum* is best known for it's role in periodontitis ([Signat et al., 2011](http://www.ncbi.nlm.nih.gov/pubmed/21220789)), but it doesn't always cause disease. Surprisingly, it is very common in the oral cavity, but rare in the gut - we found it present in 100% of saliva samples, but only 13% of stool samples from the Human Microbiome Project (HMP). For this project, I am interested in the strain-level genomic variation in *F. nucleatum* that leads to discrepancies in biogeographic colonization.

##Questions

My two basic questions and the following predictions:

1. Are there genes that are required for *F. nucleatum* to colonize the gut? If so, I will find genomic variation in *F. nucleatum* between body sites.

2. Do individuals harbor different strains of *F. nucleatum*? If so, I will find genomic variation in *F. nucleatum* between individuals.

To make the comparisons, I will need to obtain sequences of multiple *F. nucleatum* strains from the oral cavity and gut of different individuals. Here is outline of my current strategy for assembling, extracting, and comparing these genomes.

##Aquire data and pool samples

The HMP performed whole-genome shotgun sequencing on the microbiomes of 162 healthy individuals with 1253 total samples including 953 from stool and oral sites ([sequence files can be found here](http://www.hmpdacc.org/HMASM/#data)). I plan to pool together multiple samples to assemble and extract the desired genomes, but I still need to determine how many samples I will need to combine. I cannot extract a genome from a single sample, because Fuso is normally a relatively rare bug and algorithms for assembling a single genome from metagenomic data aren't very sensitive. However, the fastq files from the HMP are gigantic (~20Gb each!), so I won't be able to assemble all samples at once. When deciding how to pool the samples to create assemblies, there are many variables to consider. 

Factors that will have to be used to separate pools:

1. Body site - samples were taken at different body sites. To answer question 1, I will pool the same body sites between different individuals. 

2. Subject ID - All subjects were sampled at 18 body sites, but not all samples were sequenced and passed the quality filtering step. Because of this there are some subjects who have no oral samples and some that have as many as 18. To answer question 2, I will pool different body sites in the same individual.

3. Phase - the data was sequenced and released in two phases. I cannot pool together sequences from different phases. 

4. Study day - data was collected on different days. Some individuals are sampled on different days. These cannot be pooled together. (**but can I pool together different individuals on different days?**)

**I am going to do a preliminary assembly with 14 samples that were taken on study day 1 from the tongue of 14 different individuals. If this assembly works, I will pool the samples I have into groups of about 14 based on body site to assemble.**

##Making the data manageable

1. Download samples from a single pool into working directory

2. Quality trimming - not sure if I have to do this since the HMP went through a trimming step

3. Rename sequences to include subj id, so samples can be separated later

4. Interleave paired-end reads using khmer interleave-reads.py

5. Cat to one paired-end reads file and one singleton file

6. Digital normalization on paired-end reads and singleton reads

7. Filter-abund.py?

##Assembly

By doing an assembly of the CONCOCT mock sequence data on different assemblers (see previous posts), I will assemble the normalized pooled sequences using Megahit, IDBA, and an iterative assembly in velvet. 

##CONCOCT

See my [previous post](http://agelmore.github.io/2014/11/29/CONCOCTpipeline.html) outlining the CONCOCT pipeline that I will use to bin the contigs based on sequence composition and coverage. 


##Classify genomes 

After I separate the contigs into bins, I can use Fuso reference genomes to pull out the bins specific to Fuso using TAXAassign or BLAST. I can do another round of assembly to complete these contigs into full-length genomes.

##Genome annotation

KEGG automated annotation server (KAAS), RAST, BASys

##Compare genomes

Once I have the complete genomes annotated from each pool, I can compare across pools (between individuals and body sites). 

##Final Biogeographical Survey

I can use completed genomes from low complexity body sites (oral cavity) to detect the presence of different strains of Fuso in high complexity body sites (stool). By using this method I won't have to assemble and extract Fuso from every sample pool, but I can use the sequences that I have already assembled as reference genomes to assemble unknown pools. By doing this, I will hopefully be able to create a table including the strains present in each individual at each body site. From this I will be able to make a full strain-level biogeographical comparison.