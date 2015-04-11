---
layout: post
title:  "Velvet reference assembly"
date:   2015-4-10
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

~~~~


Now we can run velvet by inputing the reference sequence and the sorted sam file. 

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/velvet

velveth k21 21 -reference /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single.fna -short -sam velvet.fusodb.sam

~~~~

The assembly didn't work. Here's the logfile:

~~~~
[0.000000] Reading FastA file /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database
/fuso.single.fna;
[5.576443] 1 sequences found
[5.576448] Done
[5.721275] Reading FastQ file /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/ru
n2/DN/All.D1.Tongue.run2.norm.fq;
[996.303636] 232311011 sequences found
[996.303639] Done
[996.303639] Reference mapping counters
[996.303640] Name	Read mappings
[996.303640] gi|19703352|ref|NC_003454.1|	0
[996.303640] WARNING: None of your read mappings recognized the reference sequen
ce!
[996.303641] Double check that the names are identical between reference fasta h
eaders and SAM/BAM sequences.
[1003.752174] Reading read set file k21/Sequences;
[1403.315588] 232311012 sequences found
[1403.319024] Read 1 of length 32793, longer than limit 32767
[1403.319027] You should modify recompile with the LONGSEQUENCES option (cf. man
ual)
qsub working directory absolute is
/mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/velvet
~~~~







