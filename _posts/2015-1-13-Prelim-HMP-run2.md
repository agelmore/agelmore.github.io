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

{% highlight bash %}

#!/bin/bash

khmerEnv
for i in $HMP/D1.tongue/run2/wget/*.tar.bz2; do
	tar jxvf $i
	base=(basename $i .tar.bz2)
	if [ -d $base ]; then cd $base; fi
	cat $i.denovo_duplicates_marked.trimmed.1.fastq $i.denovo_duplicates_marked.trimmed.2.fastq >> $HMP/D1.tongue/run2/cat/All.D1.tongue.run2.cat.fq
	echo $base >> $HMP/D1.tongue/run2/cat/progress.txt
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






