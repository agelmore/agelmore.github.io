---
layout: post
title:  "Bowtie assembly"
date:   2015-3-30
comments: true
---

I wasn't getting good assembly of Fusobacterium with the **de novo** assembly, so I've decided to try an assembly using a reference genome. Briefly, I will use bowtie to extract reads (from the HMP dataset) that align to Fusobacterium (from a reference database). 

I want to start by re-making my Fuso database. I'm pretty sure some draft genomes have been added since I first made the database last summer. I also want to have a pipeline so I can do this quickly in the future.

~~~~

Code to make database

~~~~

Once I make the database, I'll map the raw reads to the database using bowtie.
