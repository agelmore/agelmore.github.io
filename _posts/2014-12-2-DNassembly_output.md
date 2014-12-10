---
layout: post
title:  "Assembly stats"
date:   2014-12-6
comments: true
---

Here are the assembly statistics for Velvet, Ray, and Megahit when using the interleaved read file run through digital normalization (k=21, C=20).

If you forgot, here's the command I used to do the DN:

~~~~
python2.7 normalize-by-median.py -C 20 -k 21 -p -x 1e9 All_code_interleaved.fasta -o All.code.normalized.fasta
~~~~

###Velvet

The assembly command to refresh your memory:

{% highlight bash %}

velveth velveth_k31_code 31 -shortPaired -fasta fasta/All.code.normalized.fasta
velvetg velveth_k31_code -cov_cutoff auto -read_trkg yes

{% endhighlight %}

Notice I added the command -read_trkg. I discovered that whenever I assembled with velvet I would get a summary that would say the assembly used "0/132585322 reads", even though there was a perfectly good assembly. If you use the read tracking option, velvet will tell you how many reads it actually used. 


Tail of logfile:

{% highlight bash %}

Final graph has 366968 nodes and n50 of 9396, max 181525, total 325111290, using
 127930417/132585322 reads
 
{% endhighlight %}

Assembly statistics:

~~~~
N50: 9075
N90: 526
total contigs: 254548
average length: 1303 bp
trimmed average length: 1303 bp
greater than or equal to 100:  210168
shortest conting: 61 bp
longest contig: 181555 bp
total length: 331.880719 Mb
Percent reads used: 96.5%
~~~~

###Ray

The assembly command:

{% highlight bash %}

./Ray -k 31 -i $CONCOCT_SPECIES/run1/fasta/All.code.normalized.fasta -o $CONCOCT_SPECIES/run1/ray9

{% endhighlight %}

I'm having a lot of trouble getting this to work. Sometimes it starts over, and then it says that my output folder already exists. Not sure why this is happening. Plus, I don't think there's really a benefit to using Ray since we don't have an mpi on axiom that can run on parallel nodes. I'll put this assembler on the back burner for now since I have so many others to work with.

###Megahit

The assembly command:

~~~~
python ./megahit -m 45e9 -r $CONCOCT_SPECIES/run1/fasta/All.code.normalized.fasta --cpu-only -l 100 -o $CONCOCT_SPECIES/run1/megahit_DN
~~~~

Tail of logfile:

~~~~
Total length: 135275096, N50: 69105, Mean: 11396, number of contigs: 11870
Maximum length: 727277
Done! Time elapsed(sec.): 33.450353
[Sun Nov 30 20:12:27 2014]: Merging to output final contigs..
[Sun Nov 30 20:12:36 2014]: ALL DONE.
~~~~

Took about 2.5 hours.

Stats:

~~~~
N50: 69241
N90: 9486
total contigs: 29397
average length: 11888 bp
trimmed average length: 11864 bp
greater than or equal to 100:  29397
shortest conting: 200 bp
longest contig: 727277 bp
total length: 349.496466 Mb
~~~~

###Iterative assembly

Kathryn is so nice and pretty and ran an iterative assembly using velvet for me! If you're interested in how she does this, look [here](https://kdiverson.github.io/2014/12/03/iterative-assemblies.html).

~~~~
132585322 (100.00%) were unpaired; of these:
4790404 (3.61%) aligned 0 times
98633353 (74.39%) aligned exactly 1 time
29161565 (21.99%) aligned >1 times
N50: 10038
N90: 100
total contigs: 2097980
average length: 282 bp
trimmed average length: 281 bp
greater than or equal to 100:  1437022
shortest conting: 41 bp
longest contig: 1617347 bp
total length: 592.706371 Mb
~~~~

###Summary statistics 

How to calculate the stats:

I used a script from Kathy (located here /mnt/EXT/Schloss-data/bin/contigStats.py and to calc N50 here /share/scratch/bin/calcN50N90.pl), but you can get the stats individually using these one-liners:

~~~~
grep -c '>' final.contigs.fa #find number of contigs

awk '{/>/&&++a||b+=length()}END{print b/a}' final.contigs.fa #find average sequence length

awk '!/^>/ {next} {getline s} length(s) >= 1000 { print $0 "\n" s }' contigs.fa > contigs.1000.fa; grep -c '>' contigs.1000.fa #find number of contigs > 1kb #find number of contigs greater than 1kb and save to a new file

~~~~



Assembler | Number of contigs | N50 | N90 | Average length | Contigs > 1kb | percent of reads used | contig file name
:--------|:--------:|:--------:|:--------:|:------------:|--------:
Velvet | 254548 | 9075 | 526 | 1303 | 0 |    96.5% | velveth_k31_code/contigs.fa
Megahit | 29397 | 69241 | 9486 | 11888 |    15167 | x% | megahit_DN/final.contigs.fa
Iterative assembly | 2097980 | 10038 | 100 | 282 | ? | 96.4% 



