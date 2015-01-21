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

**In retrospect, I don't want to delete the separated read files because I will want them to map to the contigs later.**

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

File info for all the lumped samples from Day 1 tongue samples just so I can remember what they are. I'd like to keep this naming scheme pretty standard for my pipeline.

File name | directory | paired? | normalized? | size | number of reads 
:---------------|:--------:|:--------:|:--------:|:--------:|:------------:
All.pe.2.fq | cat | yes | no | 217G |  1006343072
All.D1.Tongue.sin.fq | cat | singles | no | 18G | 76037114
All.D1.Tongue.normalized.cat.fq | normalize | both | yes | 50G | 234077640
All.D1.Tongue.normalized.pe.fa | normalize | yes | yes | 24G | 188980571
All.D1.Tongue.normalized.cat.fa | normalize | both | yes | 29G | 234077640
All.code.normalized.fasta | from CONCOCT mock | yes | yes | 15G | 




** number of lines in a fastq file: cat file.fastq | echo $((`wc -l`/4))


##Megahit

Megahit doesn't use paired-end read information, so I can use the file with lumped together fastq reads.

~~~~
python ./megahit -m 45e9 -r $HMP/D1.tongue/fasta/normalize/All.D1.Tongue.normalized.cat.fq --cpu-only -l 100 -o $HMP/D1.tongue/fasta/megahit
~~~~

This assembly finished! And the assembly statistics (see [previous post](http://agelmore.github.io/2014/12/06/DNassembly_output.html) to see how I did this:

~~~~
31880485 (13.62%) aligned 0 times
111646306 (47.70%) aligned exactly 1 time
90550849 (38.68%) aligned >1 times
86.38% overall alignment rate
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

This command exceeded the MEM usage hard limit. I tried again starting the iteration at a higher k value. The smaller the k, the larger the graph has to be.

~~~~
$IDBA/idba -r $HMP/D1.tongue/fasta/normalize/All.D1.Tongue.normalized.pe.fa -o $HMP/D1.tongue/fasta/idba --mink 30
~~~~


Still can't assemble. The job gets deleted because it goes over the MEM usage hard limit. Will come back to this later and use the megahit assembly for now.



## Assembly comparison

Assembler | kmer length | Number of contigs | N50 | N90 | Average length | Contigs > 1kb | percent of reads used | assembly file name
:---------------|:--------:|:--------:|:--------:|:--------:|:------------:|:------------:|:------------:|--------:
Megahit (non-paired) | iterative (21-99, step 2) | 1181605 | 606 | 243 | 515 |  99223 | 86.38% | megahit/final.contig.fa
IDBA | iterative (30, 40, 50) | Number of contigs | N50 | N90 | Average length | Contigs > 1kb | percent of reads used | assembly file name


##CONCOCT pipeline

Adapted from [previous post](http://agelmore.github.io/2014/11/29/CONCOCTpipeline.html).

~~~~
cd $HMP/D1.tongue/fasta
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m megahit/final.contigs.fa > megahit/megahit.contigs_c10K.fa
~~~~

Before I can map reads to the contigs, I have to lump the pe and sin files together for each sample in the cat folder. To do that I made this quick one-liner for loop:

~~~~
for i in SRS013502 SRS013705 SRS013818 SRS014573 SRS015174 SRS015434 SRS015537 SRS016342 SRS016501 SRS016529 SRS016569 SRS017209; do zcat $i* > $i.cat.fq; done

gzip *.cat.fq
~~~~

Okay now map the reads to the contigs:
final.contigs.fa.bowtie
~~~~
#!/bin/bash
for f in $HMP/D1.tongue/fasta/cat/Sample*R1*.fasta.gz; do 
  	gunzip -c $f > $(echo $f | sed s/".gz"/""/)
	gunzip -c $(echo $f | sed s/"R1"/"R2"/) > $(echo $(echo $f | sed s/".gz"/""/) | sed s/"R1"/"R2"/)
	mkdir -p $CONCOCT_SPECIES/map/$(basename $f .gz)
	cd $CONCOCT_SPECIES/map/$(basename $f .gz)
	$CONCOCT/scripts/map-bowtie2-markduplicates.sh -ct 1 -p '-f' $(echo $f | sed s/".gz"/""/) $(echo $(echo $f | sed s/".gz"/""/) | sed s/"R1"/"R2"/) pair $CONCOCT_SPECIES/run1/assembler/assembler.contigs_c10K.fa asm_assembler bowtie2_assembler
	rm $(echo $f | sed s/".gz"/""/)
	rm $(echo $(echo $f | sed s/".gz"/""/) | sed s/"R1"/"R2"/)
	cd ../..
done
~~~~




