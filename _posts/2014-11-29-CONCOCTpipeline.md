---
layout: post
title:  "CONCOCT pipeline"
date:   2014-11-29
comments: true
---

CONCOCT was recently published in *Nature methods* by [Alneburg et al.](http://www-ncbi-nlm-nih-gov.proxy.lib.umich.edu/pubmed/?term=binning+metagenomic+contigs+by+coverage+and+composition) They introduce a new metagenomic assembly method that bins contigs based on read coverage and kmer frequency information. Here is a basic pipeline that I've put together with notes from the [read the docs](https://concoct.readthedocs.org/en/latest/).

This method uses two pieces of information to bin contigs - coverage and composition. Coverage is found by mapping reads onto the contigs.  Composition is found by calculating the frequency of each kmer (usually k=4) in each contig. The coverage and composition vectors are combined and the contigs clustered based on the combined vector.

##Outline

1.	Assemble reads into contigs (Velvet, Ray, Megahit, SPAdes, IDBA?)
2.	Process contigs (Cut into 10kb chunks, remove duplicates)
3.	Map reads to contigs
4.	Create coverage table
5.	Run concoct to bin contigs
6.	Package includes scripts to generate plots

##1. Assembly

Working on this in other posts.

##2. Process contigs

First, cut up the contigs into chunks less than 10kb so that there isn't a bias towards mapping onto long contigs. 

**To run concoct on axiom, I have to enter the concoct environment that I created during installation. I made an alias so I just type "concoctenv" before submitting my jobs.**

```
cd $CONCOCT_SPECIES/run1
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m assembler/final.contigs.fa > assembler/assembler.congits_c10K.fa
```

Then, use bowtie2-build to make these contigs into a reference for the reads to be mapped.

```
bowtie2-build assembler/assembler.congits_c10K.fa assembler/assembler.congits._c10K.fa
```

##3. Map reads to contigs

The next step is to map all of the reads to the assembled contigs. I was having a hard time trouble shooting this huge for loop, so I put the whole thing into a bash script (/Users/Amanda/Documents/Schloss/Fuso/concoct/paperdata/bash script/markdup_bash4.sh) so I could run quicksubmit "sh markdup_bash4.sh" to submit the job. The script will unzip the reads one by one and then run the CONCOCT map-bowtie2-markduplicates.sh script. In order to run this script you must have bowtie2, BEDtools, and piccard tools installed. It requires that you create a variable in your bash_profile that points to the MarkDuplicates piccard tool.

```
export MRKDUP=/home/username/src/picard-tools-1.77/MarkDuplicates.jar
```

{% highlight bash %}

#!/bin/bash
for f in $CONCOCT_SPECIES/Sample*R1*.fasta.gz; do 
  	gunzip -c $f > $(echo $f | sed s/".gz"/""/)
	gunzip -c $(echo $f | sed s/"R1"/"R2"/) > $(echo $(echo $f | sed s/".gz"/""/) | sed s/"R1"/"R2"/)
	mkdir -p $CONCOCT_SPECIES/map/$(basename $f .gz)
	cd $CONCOCT_SPECIES/map/$(basename $f .gz)
	$CONCOCT/scripts/map-bowtie2-markduplicates.sh -ct 1 -p '-f' $(echo $f | sed s/".gz"/""/) $(echo $(echo $f | sed s/".gz"/""/) | sed s/"R1"/"R2"/) pair $CONCOCT_SPECIES/run1/assembler/assembler.congits_c10K.fa asm_assembler bowtie2_assembler
	rm $(echo $f | sed s/".gz"/""/)
	rm $(echo $(echo $f | sed s/".gz"/""/) | sed s/"R1"/"R2"/)
	cd ../..
done

{% endhighlight %}

##4. Create coverage table

Once the BAM tables have been made, use the gen_input_table.py script to make a single table the includes all the coverage information from each sample per contig. When the previous script mapped the reads, it created a separate BAM file for each sample in the folder called map. This script takes all of that coverage information and cats it to a single tsv file. Put that file into a folder called concoct-input where all of the concoct input scripts will go. 

{% highlight bash %}

cd $CONCOCT_SPECIES/map
python $CONCOCT/scripts/gen_input_table.py --isbedfiles --samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) $CONCOCT_SPECIES/run1/assembler/assembler.congits_c10K.fa */bowtie2_assembler/asm_assembler_pair-smds.coverage> concoct_inputtable.tsv
mkdir $CONCOCT_SPECIES/run1/concoct-input
mv concoct_inputtable.tsv $CONCOCT_SPECIES/run1/concoct-input/

{% endhighlight %}

Then cut the file to remove unneeded columns. 

```
cut -f1,3-26 concoct-input/concoct_inputtable.tsv > concoct-input/concoct_inputtableR.tsv
```

##5. Run concoct to bin contigs

Finally, run concoct! The -c parameter tells the clustering algorithm where to start, but isn't the absolute number of bins that will be made. However, this is a parameter that we could go back and play with later. All of the output gets put into the concoct-output file including a logfile that tells you how the run went

```
concoct -c 40 --coverage_file concoct-input/concoct_inputtableR.tsv --composition_file $CONCOCT_SPECIES/run1/assembler/assembler.congits_c10K.fa -b concoct-output/
```
