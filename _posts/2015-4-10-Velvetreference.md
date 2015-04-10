---
layout: post
title:  "Velvet reference assembly"
date:   2015-3-30
comments: true
---
The assembler [Velvet](http://www.ebi.ac.uk/~zerbino/velvet/Columbus_manual.pdf) has the option of including a reference sequence to aid in the assembly. Thought I might as well give it a try. 

The reference genome must be a single genome with non overlapping regions. I picked the fully assembled genome `gi|19703352|ref|NC_003454.1| Fusobacterium nucleatum subsp. nucleatum ATCC 25586 chromosome, complete genome` from NCBI.

~~~~

cd /mnt/EXT/Schloss-data/amanda/Fuso/extract/Fusobacterium_nucleatum_ATCC_25586_uid57885

cp * ../Database
cd ../Database
mv NC_003454.fna fuso.single.fna
bowtie2-build fuso.single.fna fuso.single
~~~~


##Map reads

First step is the use bowtie to map the reads to the reference and generate a SAM file. I'll use the normalized reads from the tongue samples (`/mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/run2/DN/All.D1.Tongue.run2.norm.fq`), because velvet won't be able to assemble with the full data set.

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/velvet

bowtie2 /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single -q /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/run2/DN/All.D1.Tongue.run2.norm.fq -p 16 -S velvet.fusodb.sam  

#sort the alignment
sort velvet.fusodb.sam > velvet.fusodb.sorted.sam
~~~~


Now we can run velvet by inputing the reference sequence and the sorted sam file. 

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/velvet

velveth k21 21 -reference /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single.fna -short -sam velvet.fusodb.sorted.sam

~~~~
