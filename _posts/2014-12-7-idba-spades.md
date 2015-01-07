---
layout: post
title:  "IDBA and SPAdes"
date:   2015-1-4
comments: true
---

I decided to try assembling with IDBA and SPAdes and comparing them to my velvet and megahit assemblies. Updated with complete assembly statistics.

###IDBA

IDBA is another de Bruijn graph assembler that uses short reads. What makes it better than Velvet is that it does an iterative assembly starting at a small k and building with larger k values. This requires more computational power, but should create a better assembly. I downloaded IDBA from [here](http://i.cs.hku.hk/~alse/hkubrg/projects/idba_ud/index.html). There is also an option to use a version of IDBA that is good with uneven depth metagenomic reads (IDBA-UD), but since I will be assembling using data that has been digitally normalized, I won't use this version.

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/idba
wget http://hku-idba.googlecode.com/files/idba-1.1.1.tar.gz
tar -zxf idba-1.1.1.tar.gz
cd idba-1.1.1
./configure
make
~~~~

To run the assembly:

~~~~
$IDBA/idba -r fasta/All.code.normalized.fasta -o idba
~~~~

It finished and created a contig.fa and a scaffold.fa file. The second uses linkage info from the paired-end reads. Here are the stats on the scaffolds from the logfile as well as the stats from Kathryn's script. Nice to see that they're the same.

~~~~
contigs: 121655 n50: 73061 max: 1064400 mean: 2899 total length: 352768328 n80: 
23533

1506871 (1.14%) aligned 0 times
128796009 (97.14%) aligned exactly 1 time
2282442 (1.72%) aligned >1 times
98.86% overall alignment rate
N50: 73061
N90: 8658
total contigs: 121655
average length: 2899 bp
trimmed average length: 2891 bp
greater than or equal to 100:  60948
shortest conting: 50 bp
longest contig: 1064400 bp
total length: 352.768328 Mb
~~~~

###SPAdes

[SPAdes](http://www-ncbi-nlm-nih-gov.proxy.lib.umich.edu/pubmed/22506599) is another assembler that uses de Bruijn graphs and can incorporate paired read information.

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/spades
wget http://spades.bioinf.spbau.ru/release3.1.1/SPAdes-3.1.1-Linux.tar.gz
tar -zxf SPAdes-3.1.1-Linux.tar.gz
cd SPAdes-3.1.1-Linux/bin
~~~~

That's it. No installation required.

Something cool: SPAdes has a --continue option that will pick up the assembly from the most recent checkpoint if it crashes for some reason while assembling.

The assembly command:

~~~~
python $SPADES/spades.py --12 $CONCOCT_SPECIES/run1/fasta/All.code.normalized.fasta -o $CONCOCT_SPECIES/run1/spades2 --only-assembler --continue
~~~~



###Updated stats

Here is the same table including the CONCOCT paper data using Ray, idba and SPAdes assembly data:

Assembler | kmer length | Number of contigs | N50 | N90 | Average length | Contigs > 1kb | percent of reads used | assembly file name
:---------------|:--------:|:--------:|:--------:|:--------:|:------------:|:------------:|:------------:|--------:
CONCOCT paper (Ray) | 41 | 13068 | 153630 | not shown | not shown |  37627 | 99% | 
Velvet | 31 | 254548 | 9075 | 526 | 1303 | 0 |    91.6% | velveth_k31_code/contigs.fa
Megahit | iterative (21-99, step 2)| 29397 | 69241 | 9486 | 11888 |    15167 | 99.84% | megahit_DN/final.contigs.fa
Kathryn's Iterative assembly | iterative (49-21) | 2097980 | 10038 | 100 | 282 | 19079 | 96.4% 
IDBA | iterative (20, 30, 40, 50) | 121655 | 73061 | 8658 | 2899 |    12786 | 98.86% | idba/scaffold.fa
SPAdes | iterative (21,33,55) | Number of contigs | N50 | N90 | Average length | Contigs > 1kb | percent of reads used | assembly file name

