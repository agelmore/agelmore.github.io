---
layout: post
title:  "Run2 - Preliminary HMP data assembly"
date:   2015-1-22
comments: true
---

As I was running the first preliminary assembly, I came across a couple of problems with which files I was using. I'm going to start over with a few more samples and fix those problems with this pipeline.

First of all, what to change:

1. Don't use the singletons file in the megahit assembly. Megahit can't use paired-end read information, but it makes things easier later if I ignore these.

2. Don't need to do DN on the singletons.

3. Don't delete the separated read files because I need them to do CONCOCT.

4. Use more samples (20+)


###Download

I used a few more samples this time so that I had a total of 20. They are still all from different individuals from the tongue and collected on the first day for that individual. Here are the locations for wget.

{% gist 96674047b7c23dfdf92c %}

###Managing the data

~~~~
tar zxvf SRS*
~~~~


{% highlight bash %}

#!/bin/bash

for i in $HMP/D1.tongue/run2/wget/SRS*; do
	cd $i
	cat *.denovo_duplicates_marked.trimmed.1.fastq *.denovo_duplicates_marked.trimmed.2.fastq >> $HMP/D1.tongue/run2/cat/All.D1.tongue.run2.cat.fq
	echo $i >> $HMP/D1.tongue/run2/cat/progress.txt
	gzip *q
	cd $HMP/D1.tongue/run2/wget
done

{% endhighlight %}

###Digital normalization

~~~~
khmerEnv
cd $HMP/D1.tongue/run2
mkdir DN
cd DN
python2.7 /mnt/EXT/Schloss-data/amanda/Fuso/khmer/khmerEnv/bin/normalize-by-median.py -C 20 -k 21 -x 1e9 ../cat/All.D1.tongue.run2.cat.fq -s All.D1.Tongue.run2.savetable -o All.D1.Tongue.run2.norm.fq
~~~~

###Megahit

~~~~
cd $HMP/D1.tongue/run2
mkdir megahit
cp /mnt/EXT/Schloss-data/amanda/Fuso/megahit/megahit ../../HMP/D1.tongue/run2/megahit/
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/run2/megahit/
python ./megahit -m 45e9 -r $HMP/D1.tongue/run2/DN/All.D1.Tongue.run2.norm.fq --cpu-only -l 100 -o $HMP/D1.tongue/run2/megahit
~~~~

###Assembly stats

~~~~
python /mnt/EXT/Schloss-data/bin/contigStats.py $HMP/D1.tongue/run2/megahit/final.contigs.fa

perl /share/scratch/bin/calcN50N90.pl $HMP/D1.tongue/run2/megahit/final.contigs.fa

bowtie2-build $HMP/D1.tongue/run2/megahit/final.contigs.fa $HMP/D1.tongue/run2/megahit/final.contigs.fa.bowtie
bowtie2 final.contigs.fa.bowtie -q ../DN/All.D1.Tongue.run2.norm.fq -p 16 -S megahit.aligned.sam 


cd $HMP/D1.tongue/run2/megahit/; awk '!/^>/ {next} {getline s} length(s) >= 1000 { print $0 "\n" s }' final.contigs.fa > final.contigs.1000.fa; grep -c '>' final.contigs.1000.fa 
~~~~

Output:

~~~~
35870276 (15.44%) aligned 0 times
    113973643 (49.06%) aligned exactly 1 time
    82467092 (35.50%) aligned >1 times
84.56% overall alignment rate
N50: 587
N90: 24
total contigs: 1403622
average length: 507 bp
trimmed average length: 507 bp
greater than or equal to 100:  1403622
shortest conting: 200 bp
longest contig: 259539 bp
total length: 712.681025 Mb
contigs > 1kb: 111168
~~~~



Assembler | kmer length | Number of contigs | N50 | N90 | Average length | Contigs > 1kb | percent of reads used | assembly file name
:---------------|:--------:|:--------:|:--------:|:--------:|:------------:|:------------:|:------------:|--------:
Megahit (non-paired) | iterative (21-99, step 2) | 1403622 | 587 | 24 | 507 |  111168 | 84.56% | $HMP/D1.tongue/run2/megahit/final.contigs.fa

###CONCOCT pipeline

~~~~
concoctenv
cd $HMP/D1.tongue/run2
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m megahit/final.contigs.fa > megahit/megahit.contigs_c10K.fa

bowtie2-build megahit/megahit.contigs_c10K.fa megahit/megahit.contigs_c10K_bowtie.fa
~~~~

Map reads to contigs. Saved in bash script called map.reads.sh

{% highlight bash %}

#!/bin/bash

for f in $HMP/D1.tongue/run2/wget/SRS*; do 
  	cd $f;
  	gunzip *;
  	mkdir -p $HMP/D1.tongue/run2/concoct/map/$(basename $f);
    cd $HMP/D1.tongue/run2/concoct/map/$(basename $f);
    $CONCOCT/scripts/map-bowtie2-markduplicates.sh -ct 1 -p '-q' $f/$(basename $f).denovo_duplicates_marked.trimmed.1.fastq $f/$(basename $f).denovo_duplicates_marked.trimmed.2.fastq pair $HMP/D1.tongue/run2/megahit/megahit.contigs_c10K.fa asm bowtie2;
    gzip $f/*.fastq
done


{% endhighlight %}

Okay that was taking a lot time so I submitted it again with each sample as a separate job. I'm sure I can do this later by making a file that will submit jobs in a for loop. I bet Kathryn has showed me how to do that before. For now I created individual sh files and submitted them all manually.

{% highlight bash %}

#!/bin/bash

f=$HMP/D1.tongue/run2/wget/SRS
cd $f;
gunzip *;
mkdir -p $HMP/D1.tongue/run2/concoct/map/$(basename $f);cd $HMP/D1.tongue/run2/concoct/map/$(basename $f);
$CONCOCT/scripts/map-bowtie2-markduplicates.sh -ct 1 -p '-q' $f/$(basename $f).denovo_duplicates_marked.trimmed.1.fastq $f/$(basename $f).denovo_duplicates_marked.trimmed.2.fastq pair $HMP/D1.tongue/run2/megahit/megahit.contigs_c10K.fa asm bowtie2;
gzip $f/*.fastq


{% endhighlight %}

For some reason the script didn't work for all of the samples. A few of them didn't make the smds.coverage files, but did make the sam files. I removed those samples to figure that out later, and went on to make the coverage table with just 13 samples. 

Generate coverage table:

~~~~
cd $HMP/D1.tongue/run2/concoct/map
python $CONCOCT/scripts/gen_input_table.py --isbedfiles --samplenames <(for s in SRS*; do echo $s; done) ../../megahit/megahit.contigs_c10K.fa SRS*/bowtie2/asm_pair-smds.coverage > concoct_inputtable.tsv
mkdir $HMP/D1.tongue/run2/concoct/concoct-input
mv concoct_inputtable.tsv $HMP/D1.tongue/run2/concoct/concoct-input

~~~~

Generate linkage table:

~~~~
cd $HMP/D1.tongue/run2/concoct/map
python $CONCOCT/scripts/bam_to_linkage.py -m 8 --regionlength 500 --fullsearch --samplenames <(for s in SRS*; do echo $s; done) ../../megahit/megahit.contigs_c10K.fa SRS*/bowtie2/asm_pair-smds.bam > concoct_linkage.tsv
mv concoct_linkage.tsv $HMP/D1.tongue/run2/concoct/concoct-input

~~~~

Cut table:

~~~~
cd $HMP/D1.tongue/run2/concoct
cut -f1,3-26 concoct-input/concoct_inputtable.tsv > concoct-input/concoct_inputtableR.tsv
~~~~

Run CONCOCT:

~~~~
concoct -c 40 --coverage_file concoct-input/concoct_inputtableR.tsv --composition_file $HMP/D1.tongue/run2/megahit/megahit.contigs_c10K.fa -b concoct-output/
~~~~

Concoct finished in less than a day. I ran the ClusterPlot.R script that comes with the CONCOCT software and generated a scatterplot. It doesn't look like much.


I combined the file with the cluster numbers (clustering_gt1000.csv) with the blast output (blast.fuso.all) which should be able to tell me where the contigs that had fuso clustered to. This file only has 3580 contigs, probably because I didn't use all the samples that I had done the assembly with when I ran concoct. Looking for unique clusters in this file (`cat blast.concoct.merged.txt | cut -f 2 | sort | uniq`) showed that almost every cluster was found in these contigs. Something is wrong here. 



