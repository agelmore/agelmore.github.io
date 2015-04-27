---
layout: post
title:  "CONCOCT on reference assembly"
date:   2015-4-27
comments: true
---

The SCG analysis and coverage map show that I may have a couple different genomes in this assembly. The CONCOCT paper showed that they could get separation at the strain level, so I'm going to try it to see if I can get clusters of individual strains. 

It might not work with only 20 samples, but it's worth a shot.

Assembled contigs: `$HMP/D1.tongue/reference/bowtie/DN/megahit2/final.contigs.fa`

Directory to run concoct: `/mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/DN/megahit2/concoct
`

Pipeline:

**Map reads to contigs**

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/DN/megahit2/concoct
concoctenv

mkdir map
cd map
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