---
layout: post
title:  "Preliminary HMP data assembly"
date:   2015-1-12
comments: true
---

To choose the best sample size for pooling body sites and individuals, I did a preliminary assembly of 14 samples from the tongue of 14 different individuals taken on study day 1.

###Download

I downloaded the following samples to axiom using wget:

{% gist 678041d9daa5cf164bfe %}

###Managing the data

I wanted to interleave the paired-end reads files, but the files are too big to unzip all at once. I wrote this shell script to extract the files one by one, interleave read 1 and 2, move and zip the interleaved and singleton files to another directory, and then delete the unzipped directory.

{% highlight bash %}

#!/bin/bash

for i in $HMP/D1.tongue/fasta/*.tar.bz2; do
	tar jxvf $i
	cd $(basename $i .tar.bz2)
	python2.7 /mnt/EXT/Schloss-data/amanda/Fuso/khmer/khmerEnv/bin/interleave-reads.py *.1.fastq *.2.fastq -o $(basename $i .tar.bz2).pe.fq
	mv $(basename $i .tar.bz2).pe.fq ../cat
	mv *.singleton.fastq ../cat
	cd ..
	rm -r $(basename $i .tar.bz2)
	cd cat
	gzip *q
	cd ..
done

{% endhighlight %}

It finished in a day. I then combined all of the pe read files into a single file that could be used for assembly.

~~~~

zcat SRS*.pe.fq > All.D1.Tongue.pe.fq
zcat *.singleton.fastq > All.D1.Tongue.sin.fq

~~~~

I think it finished but how can I tell if it's complete? Count number of reads in each of the files and then count the reads in the cat file? Or just do the cat again.

###Digital Normalization

I'm waiting on the files to merge, but I plan to use these parameters to run digital normalization when it finishes.

~~~~
python2.7 /mnt/EXT/Schloss-data/amanda/Fuso/khmer/khmerEnv/bin/normalize-by-median.py -C 20 -k 21 -x 1e9 ../cat/All.pe.2.fq -s All.D1.Tongue.pe.savetable -o All.D1.Tongue.normalized.pe.fq
~~~~

Output:

~~~~
DONE with ../cat/All.D1.Tongue.sin.fq; kept 45097069 of 76037114 or 59%
output in All.D1.Tongue.normalized.sin.fq
Saving k-mer counting table through ../cat/All.D1.Tongue.sin.fq
...saving to All.D1.Tongue.sin.savetable
fp rate estimated to be 0.113

DONE with ../cat/All.pe.2.fq; kept 188980571 of 1006343072 or 18%
output in All.D1.Tongue.normalized.pe.fq
Saving k-mer counting table through ../cat/All.pe.2.fq
...saving to All.D1.Tongue.pe.savetable
fp rate estimated to be 0.757
~~~~

##Assemblies

I'm not sure what to do with the paired-end reads and singles, so I started by lumping them together and assembling like they're all unpaired.

##Megahit

~~~~
python ./megahit -m 45e9 -r $HMP/D1.tongue/fasta/normalize/All.D1.Tongue.normalized.cat.fq --cpu-only -l 100 -o $HMP/D1.tongue/fasta/megahit
~~~~

Can I assemble with paired reads? If so, what should I do with the singles?

This assembly finished! And the assembly statistics (see [previous post](http://agelmore.github.io/2014/12/06/DNassembly_output.html) to see how I did this:

bowtie2-build final.contigs.fa final.contigs.fa.bowtie
bowtie2 final.contigs.fa.bowtie -q ../normalize/All.D1.Tongue.normalized.cat.fq -p 16 -S megahit.1.aligned.sam

~~~~
N50: 606
N90: 243
total contigs: 1181605
average length: 515 bp
trimmed average length: 515 bp
greater than or equal to 100:  1181605
shortest conting: 200 bp
longest contig: 201960 bp
total length: 609.704148 Mb
contigs > 1kb: 99223
~~~~


##IDBA

IDBA only works with paired-end reads, so I had to omit all the singletons. I wonder if there is a way to assemble with them as well, but I can't with the package I have. I also had to use the fq2fa script to change the file into a fasta file. 

~~~~
$IDBA/idba -r $HMP/D1.tongue/fasta/normalize/All.D1.Tongue.pe.fa -o $HMP/D1.tongue/fasta/idba

~~~~

## Assembly comparison

Assembler | kmer length | Number of contigs | N50 | N90 | Average length | Contigs > 1kb | percent of reads used | assembly file name
:---------------|:--------:|:--------:|:--------:|:--------:|:------------:|:------------:|:------------:|--------:
Megahit (non-paired) | iterative (21-99, step 2) | 1181605 | 606 | 243 | 515 |  99223 | 99% | megahit/final.contig.fa
IDBA | iterative (20, 30, 40, 50) | Number of contigs | N50 | N90 | Average length | Contigs > 1kb | percent of reads used | assembly file name




