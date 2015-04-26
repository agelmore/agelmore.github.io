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

#run genomeCoverageBed. The d option produces the coverage at each base. Alternatively, the bg option will produce a BEDGRAPH output file where the coverage of each range is listed. 
bedtools genomecov -d -ibam bowtie.single.sorted.bam.bam -g /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single.fna.fai > bowtie.single.density

~~~~

It finished. I pulled the file to my computer to make a plot in R.

~~~~
setwd("~/Documents/Schloss/Fuso/Assemblies/D1 tongue/bedgraph")
x<- read.table(file='bowtie.single.density')
plot(x$V2,x$V3, main="Coverage of bowtie extracted reads on F. nucleatum", xlab="location in genome", ylab='base coverage', type='l')

~~~~

Here is the coverage map of the reads from the tongue samples that aligned to the Fuso database. **I did not run these reads through digital normalization or an assembler**. 

![Coverage map reads extracted from single]({{ site.url }}/images/bowtie.single.coveragemap.png)


Looks like pretty even coverage except for those two spots. Next I want to see what happens to the coverage after digital normalization (which I pretty much have to use in order to assemble with big files). Also, should figure out what those high coverage spots are.

#Coverage map of assembled contigs on single F. nucleatum

I will use the megahit assembly of normalized extracted reads from the full fuso db. Because these are extracted from the full database, it's very possible that there will be contigs that don't map to the F. nucleatum genome. However, I'm going to see what the coverage looks like for just nucelatum. Also, I'll make a before and after assembly coverage map to see what the change is.


###Coverage after extraction

I have already extracted from the full Fuso database, but I want to visualize the coverage against a single genome (F. nucleatum). So I need to take my extracted fastq file and align it to the single genome, then make the coverage map. I can pipe this so I don't have to submit multiple commands.

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie

bowtie2 /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single -q /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/bowtie.fusodb.mapped.cut.fastq -p 16 -S bowtie.coverage.beforeDN.sam; samtools view -Sb bowtie.coverage.beforeDN.sam > bowtie.coverage.beforeDN.bam; samtools sort bowtie.coverage.beforeDN.bam bowtie.coverage.beforeDN.sorted; bedtools genomecov -d -ibam bowtie.coverage.beforeDN.sorted.bam -g /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single.fna.fai > maps/bowtie.coverage.beforeDN.sorted.density 
~~~~

Then in R:

~~~~
setwd("/mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/maps")
png('coverage.beforeDN.png')
x<- read.table(file='bowtie.coverage.beforeDN.sorted.density')
plot(x$V2,x$V3, main="Coverage of bowtie extracted reads on F. nucleatum", xlab="location in genome", ylab='base coverage', type='l', file='coverage.beforeDN.png')
dev.off()
~~~~

Pull to local blog:

~~~~
cd /Users/Amanda/Documents/Schloss/agelmore.github.io/images
sftp agelmore@axiom.ccmb.med.umich.edu
get /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/maps/coverage.beforeDN.png
~~~~

![Coverage map reads extracted from full database]({{ site.url }}/images/coverage.beforeDN.png)

It makes sense that this graph looks really similar to the one above. The only difference is that the reads here were extracted from the whole fuso database and then I graphed the ones that aligned to this single genome. Before they were extracted from the single genome. 

###Coverage after digital normalization

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/DN

bowtie2 /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single -q /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/DN/bowtie.fusodb.mapped.cut.normalized.fastq -p 16 -S bowtie.coverage.afterDN.sam; samtools view -Sb bowtie.coverage.afterDN.sam > bowtie.coverage.afterDN.bam; samtools sort bowtie.coverage.afterDN.bam bowtie.coverage.afterDN.sorted; bedtools genomecov -d -ibam bowtie.coverage.afterDN.sorted.bam -g /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single.fna.fai > ../maps/bowtie.coverage.afterDN.sorted.density 

setwd("/mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/maps")
png('coverage.afterDN.png')
x<- read.table(file='bowtie.coverage.afterDN.sorted.density')
plot(x$V2,x$V3, main="Coverage of bowtie extracted and normalized reads on F. nucleatum", xlab="location in genome", ylab='base coverage', type='l', file='coverage.afterDN.png')
dev.off()
~~~~

![Coverage map reads extracted from full database and normalized]({{ site.url }}/images/coverage.afterDN.png)

Cool! You can really see that the normalization works. The coverage decreases to about 200 per basepair and it seems to even out a bit. Also interesting that those two spots around .75Mbp and 1.1Mbp - **I need to check what those are**.





