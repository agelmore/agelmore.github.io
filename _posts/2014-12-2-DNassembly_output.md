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
 0/132585322 reads
 
{% endhighlight %}

###Ray

The assembly command:

{% highlight bash %}

./Ray -k 31 -i $CONCOCT_SPECIES/run1/fasta/All.code.normalized.fasta -o $CONCOCT_SPECIES/run1/ray9

{% endhighlight %}

I'm having a lot of trouble getting this to work. Sometimes it starts over, and then it says that my output folder already exists. Not sure why this is happening. Plus, I don't think there's really a benefit to using Ray since we don't have an mpi on axiom that can run on parallel nodes. I'll put this assembler on the back burner for now since I have so many others to work with.

###Megahit

The assembly command:

Tail of logfile:

~~~~
Total length: 135275096, N50: 69105, Mean: 11396, number of contigs: 11870
Maximum length: 727277
Done! Time elapsed(sec.): 33.450353
[Sun Nov 30 20:12:27 2014]: Merging to output final contigs..
[Sun Nov 30 20:12:36 2014]: ALL DONE.
~~~~

###Summary statistics 

How to calculate some stats:


~~~~
grep -c '>' final.contigs.fa #find number of contigs

awk '{/>/&&++a||b+=length()}END{print b/a}' final.contigs.fa #find average sequence length

awk '!/^>/ {next} {getline s} length(s) >= 1000 { print $0 "\n" s }' contigs.fa > contigs.1000.fa; grep -c '>' contigs.1000.fa #find number of contigs > 1kb #find number of contigs greater than 1kb and save to a new file
~~~~

To calculate [N50](https://github.com/kdiverson/seqTools/blob/master/calcN50.py) I stole a script from Kathy Iverson. 


Assembler | Number of contigs | N50 | Average length | Contigs > 1kb | contig file name
-------- | -------- | -------- | -------- | --------
Velvet | 254548 | 9396 | 1303.8 | 0 | velveth_k31_code/contigs.fa
Megahit | 11870 | 69241 | 11396 | 15167 | megahit_DN/final.contigs.fa


Uh oh, 0 contigs greater than 1kb for velvet? That doesn't seem good. Let's see what the assemblers IDBA and SPAdes can do.


