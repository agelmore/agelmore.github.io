---
layout: post
title:  "Bowtie/Blast assembly"
date:   2015-3-30
comments: true
---

I wasn't getting good assembly of Fusobacterium with the *de novo* assembly, so I've decided to try an assembly using a reference genome. Briefly, I will use bowtie to extract reads (from the HMP dataset) that align to Fusobacterium (from a reference database). 

I want to start by re-making my Fuso database. I'm pretty sure some draft genomes have been added since I first made the database last summer. I also want to have a pipeline so I can do this quickly in the future. All genomes are stored on the [ncbi website](ftp://ftp.ncbi.nih.gov/genomes/Bacteria).

###First make database with full length Fuso genomes

I downloaded all of the full length genomes from ncbi since they don't take up too much space. I combined all the fuso genomes (and full length plasmids) into a single file. Here are the sequence names:

{% gist 7cdeb4fc470cb793259c %}

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/extract

wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/*    #download all genomes
for f in Fuso*; do cat $f/*.fna >> Database/fusodb.fa; done

~~~~

###Database with draft genomes added

This is a little tricker since they're a lot of them. I made a text file that has the location of all the draft genomes on the NCBI site called FASTA_location. Each of these genomes has multiple scaffolds because they aren't complete genomes. I combined all of them into a single file. 

{% gist b8a98a9745241f8ab7e1 %}

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/extract/draft

wget -i FASTA_location
cat *.fna >> draftdb.fa
~~~~

And then combine the databases into one...

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database
cat fusodb.fa ../draft/draftdb.fa > fusodb.complete.fa
~~~~

###Index databases

Might as well index these databases for bowtie2 and blastn right now.

~~~~
bowtie2-build fusodb.complete.fa fusodb.bowtie2
makeblastdb -in fusodb.complete.fa -dbtype nucl -out fusodb.blast

~~~~

#Extract reads

Now I'm ready to extract reads that map to any of the Fuso genomes. I have all the reads from the 20 tongue samples that I've been using for practice in a file called `$HMP/D1.tongue/run2/cat/All.D1.tongue.run2.cat.fq`. I'm going to map those reads to the database using bowtie2 and blastn (We'll see which is better. I'm thinking bowtie will be). 

###Bowtie2

Let's start with bowtie2.

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie

bowtie2 /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fusodb.bowtie2 -q $HMP/D1.tongue/run2/cat/All.D1.tongue.run2.cat.fq -p 16 -S bowtie.fusodb.sam 

~~~~

Now to extract the mapped reads out of the sam file. 

~~~~
#make BAM file
samtools view -bT /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fusodb.bowtie2 bowtie.fusodb.sam > bowtie.fusodb.bam

#use the -F4 option. The -F option removes the specified FLAG. The 4 flag is unmapped reads. 
samtools view -F4 bowtie.fusodb.bam > bowtie.fusodb.mapped.bam

#cut out name, sequence, and quality from sam file
cut -f1,10,11 bowtie.fusodb.mapped.bam > bowtie.fusodb.mapped.cut.bam

awk '{print "@"$1"\n"$2"\n""\+""\n"$3}' bowtie.fusodb.mapped.cut.bam > bowtie.fusodb.mapped.cut.fastq
~~~~

Final read count after bowtie filtering: 

~~~~
cat bowtie.fusodb.mapped.cut.fastq | echo $((`wc -l`/4))

cat All.D1.tongue.run2.cat.fq | echo $((`wc -l`/4))
~~~~

`75566622` reads filtered out of `1006343072` in original file which is 7.5%. 

###Assembly

~~~~
khmerEnv
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie
mkdir DN
cd DN
python2.7 /mnt/EXT/Schloss-data/amanda/Fuso/khmer/khmerEnv/bin/normalize-by-median.py -C 20 -k 21 -x 1e9 ../bowtie.fusodb.mapped.cut.fastq -s bowtie.fusodb.mapped.cut.savetable -o bowtie.fusodb.mapped.cut.normalized.fastq

cd /mnt/EXT/Schloss-data/amanda/Fuso/megahit/megahit

python ./megahit -m 45e9 -r /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/DN/bowtie.fusodb.mapped.cut.normalized.fastq --cpu-only -l 101 -o /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/megahit

##Try the assembly without DN too

python ./megahit -m 45e9 -r /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/bowtie.fusodb.mapped.cut.fastq --cpu-only -l 101 -o /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/megahit


~~~~

Finished! Assembly stats:

~~~~
N50: 457
N90: 236
total contigs: 43760
average length: 429 bp
trimmed average length: 428 bp
greater than or equal to 100:  43760
shortest conting: 200 bp
longest contig: 8669 bp
total length: 18.780803 Mb
4865248 reads; of these:
  4865248 (100.00%) were unpaired; of these:
    289783 (5.96%) aligned 0 times
    1397484 (28.72%) aligned exactly 1 time
    3177981 (65.32%) aligned >1 times
94.04% overall alignment rate
~~~~

###Summary and direction

Pretty good assembly I'd say! The contigs aren't super long, but the total assembly is long enough for there to be multiple Fuso genomes in there. I'm going to run it through the single copy core genes analysis to see if it's a single genome or multiple. (In a [new blog post](), this one is getting long). Then, maybe I'll try running it through CONCOCT to separate out any strains/species.

###Blast

And map with Blast, too. I set the max_target_seqs parameter so that reads won't hit to more than one genome (which will probably happen a lot). Unfortunately, BLAST requires FASTA so I'll have to convert that first.

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/blast

awk '{print \">\" substr(\$0,2);getline;print;getline;getline}' All.D1.tongue.run2.cat.fq > All.D1.tongue.run2.cat.fa  #I have an alias fq2fa for this in my bash_profile

blastn -db /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fusodb.blast -query $HMP/D1.tongue/run2/cat/All.D1.tongue.run2.cat.fa -out blast.fusodb -evalue 1e-5 -outfmt 6 -num_threads 16 -max_target_seqs 1

~~~~

Finished after 3 days!

Extract the sequences out of blast output:

~~~~
cut -f1 blast.fusodb > blast.fusodb.reads
/share/scratch/schloss/mothur/mothur '#get.seqs(accnos=blast.fusodb.reads, fasta=/mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/run2/cat/All.D1.tongue.run2.cat.fa)'
mv All.D1.tongue.run2.cat.pick.fa blast.fusodb.reads.fa

~~~~

`78527266` reads filtered out of `1006343072` in original file which is 7.8%.

###Blast assembly


~~~~
khmerEnv
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/blast
mkdir DN
cd DN
python2.7 /mnt/EXT/Schloss-data/amanda/Fuso/khmer/khmerEnv/bin/normalize-by-median.py -C 20 -k 21 -x 1e9 ../blast.fusodb.reads.fa -s blast.fusodb.reads.savetable -o blast.fusodb.reads.normalized.fasta

cd /mnt/EXT/Schloss-data/amanda/Fuso/megahit/megahit

with normalize
python ./megahit -m 45e9 -r /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/blast/DN/blast.fusodb.reads.normalized.fasta --cpu-only -l 101 -o /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/blast/megahit/DN

#without normalize
python ./megahit -m 45e9 -r /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/blast/blast.fusodb.reads.fasta --cpu-only -l 101 -o /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/blast/megahit


~~~~

###Assembly stats

~~~~
total contigs: 43870
average length: 416 bp
trimmed average length: 415 bp
greater than or equal to 100:  43870
shortest conting: 200 bp
longest contig: 8537 bp
total length: 18.254674 Mb
N50: 438
N90: 234



~~~~


