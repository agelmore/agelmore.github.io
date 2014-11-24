---
layout: post
title:  "Assembling after normalization"
date:   2014-11-19
comments: true
---

Interleaving and digital normalization finished, so I started Ray and Velvet assembling with the normalized read files. These will both probably take a little while.

Output from Digital normalization:


###Ray

So now we can run Ray again with the same parameters as before. This time it will be able to assemble all the reads because the dataset is much smaller.

{% highlight bash %}

mpdboot
mpiexec -n 1 ./mnt/EXT/Schloss-data/amanda/ray/Ray-2.3.1/ray-build/Ray -k 32 -i $CONCOCT_SPECIES/run1/fasta/All.code.normalized.fasta -o $CONCOCT_SPECIES/run1/ray

{% endhighlight %}


###Velvet

[Velvet](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2952100/) is another *de novo* assembler meant for assembling short reads. The assembly is divided into two steps: velveth where the reads are prepared and graphs are made and then velvetg where the assembly is done.

{% highlight bash %}

velveth velveth_k41_code 31 -shortpaired -fasta fasta/All.code.normalized.fasta
velvetg velveth_k41_code -cov_cutoff auto

{% endhighlight %}

Stay tuned for assembly statistics.