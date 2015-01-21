---
layout: post
title:  "Assembling after normalization"
date:   2014-11-24
comments: true
---

Interleaving and digital normalization finished, so I started Ray and Velvet assembling with the normalized read files. These will both probably take a little while.

Output from Digital normalization:

{% highlight bash %}

DONE with All_code_interleaved.fasta; kept 131641458 of 1488000000 or  8%
output in All.code.normalized.fasta
fp rate estimated to be 0.250

{% endhighlight %} 

###Ray

So now we can run Ray again with the same parameters as before. This time it will be able to assemble all the reads because the dataset is much smaller. Also, now I'm running ray without the MPI.

{% highlight bash %}

./Ray -k 32 -i $CONCOCT_SPECIES/run1/fasta/All.code.normalized.fasta -o $CONCOCT_SPECIES/run1/ray8

{% endhighlight %}


###Velvet

[Velvet](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2952100/) is another *de novo* assembler meant for assembling short reads. The assembly is divided into two steps: velveth where the reads are prepared and graphs are made and then velvetg where the assembly is done.

{% highlight bash %}

velveth velveth_k31_code 31 -shortPaired -fasta fasta/All.code.normalized.fasta
velvetg velveth_k31_code -cov_cutoff auto

{% endhighlight %}

###Megahit

Might as well do the megahit assembly with the normalized reads as well. Titus says that we don't have to do DN with megahit, but the assembler didn't use all the reads last time because axiom ran out of memory. Since DN is basically removing redundant reads, we'll try and see if this works.

{% highlight bash %}

python ./megahit -m 45e9 -r $CONCOCT_SPECIES/run1/fasta/All.code.normalized.fasta --cpu-only -l 100 -o $CONCOCT_SPECIES/run1/megahit_DN

{% endhighlight %}

Stay tuned for assembly statistics.