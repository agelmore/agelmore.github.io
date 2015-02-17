---
layout: post
title:  "CONCOCT with 1kb reads"
date:   2015-2-16
comments: true
---

I reran concoct using only the contigs that are longer than 1kb, because the tetranucleatide frequencies that concoct uses to bin contigs works best with contigs greater than 1kb.

See previous post for how I filtered for the 1kb contigs.

Pipeline:

**Map reads to contigs**

~~~~
#!/bin/bash

f=$HMP/D1.tongue/run2/wget/SRS
cd $f;
gunzip *;
mkdir -p $HMP/D1.tongue/run2/concoct/map/$(basename $f);cd $HMP/D1.tongue/run2/concoct/map/$(basename $f);
$CONCOCT/scripts/map-bowtie2-markduplicates.sh -ct 1 -p '-q' $f/$(basename $f).denovo_duplicates_marked.trimmed.1.fastq $f/$(basename $f).denovo_duplicates_marked.trimmed.2.fastq pair $HMP/D1.tongue/run2/megahit/megahit.1000.contigs_c10K.fa asm bowtie2;
gzip $f/*.fastq
~~~~

**Generate input table**

~~~~
cd $HMP/D1.tongue/run2/concoct/map
python $CONCOCT/scripts/gen_input_table.py --isbedfiles --samplenames <(for s in SRS*; do echo $s; done) ../../megahit/megahit.1000.contigs_c10K.fa SRS*/bowtie2/asm_pair-smds.coverage > concoct_inputtable.tsv
mkdir $HMP/D1.tongue/run2/concoct/1kb/concoct-input
mv concoct_inputtable.tsv $HMP/D1.tongue/run2/concoct/1kb/concoct-input

~~~~

**Cut table:**

~~~~
cd $HMP/D1.tongue/run2/concoct/1kb
cut -f1,3-26 concoct-input/concoct_inputtable.tsv > concoct-input/concoct_inputtableR.tsv
~~~~

**Run CONCOCT:**

~~~~
concoct -c 40 --coverage_file concoct-input/concoct_inputtableR.tsv --composition_file $HMP/D1.tongue/run2/megahit/megahit.1000.contigs_c10K.fa -b concoct-output/
~~~~

**Run ClusterPlot.R script**

![ClusterPlot]({{ site.url }}/images/ClusterPlot1kb.png)

This plot shows how all of the contigs cluster into the bins on a PCA plot. To see how the blasted Fuso contigs clustered, I created a histogram of the frequency of fuso contigs in different CONOCCT bins. 

![Cluster Histogram]({{ site.url }}/images/ClusterPlot1kb.png)

This actually worked pretty well! 17, 21, and 25 had the most Fuso contigs. Next, I will re-assemble those 3 clusters and see if I can get some longer Fuso contigs that I can separate into strains.

