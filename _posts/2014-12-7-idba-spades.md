---
layout: post
title:  "IDBA and SPAdes"
date:   2014-12-7
comments: true
---

I decided to try assembling with IDBA and SPAdes and comparing them to my velvet and megahit assemblies.

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

I'll add the output when it finishes.

###SPAdes

[SPAdes](http://www-ncbi-nlm-nih-gov.proxy.lib.umich.edu/pubmed/22506599) is another assembler that uses de Bruijn graphs and can incorporate paired read information.

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/spades
wget http://spades.bioinf.spbau.ru/release3.1.1/SPAdes-3.1.1-Darwin.tar.gz
tar -zxf SPAdes-3.1.1-Darwin.tar.gz
cd SPAdes-3.1.1-Darwin/bin
~~~~

That's it. No installation required.

Okay, actually that didn't work because it's meant for Mac OS and not axiom, so I followed the directions to download and compile the source code:

~~~~
wget http://spades.bioinf.spbau.ru/release3.1.1/SPAdes-3.1.1.tar.gz
tar -xzf SPAdes-3.1.1.tar.gz
cd SPAdes-3.1.1
./spades_compile.sh
~~~~

The assembly command:

~~~~
python $SPADES/spades.py --12 fasta/All.code.normalized.fasta -o spades --only-assembler

~~~~

I'll add the output when it finishes.
