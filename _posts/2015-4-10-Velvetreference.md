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
[9.401484] 1 sequences found
[9.401492] Done
[9.402172] Reading SAM file velvet.fusodb.sam
[2844.523330] 232311011 reads found.
[2844.523333] Done
[2844.523335] Reference mapping counters
[2844.523337] Name	Read mappings
[2844.523337] gi|19703352|ref|NC_003454.1|	900570
[2844.643837] Reading read set file k21/Sequences;
[4138.882268] 232311012 sequences found
[4138.884069] Read 1 of length 32793, longer than limit 32767
[4138.884074] You should modify recompile with the LONGSEQUENCES option (cf. man
ual)
qsub working directory absolute is
/mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/velvet
~~~~

Looks like the problem is that the reference genome is vastly too long. It's 2174500 bp and it looks like the limit is 32767bp. It says I can recompile velvet with the LONGSEQUENCES option. Maybe this function is really meant to assemble certain genes or parts of the genome and not the whole thing. 

###Gene reference

Because velvet can't handle the really long reference, I hoped to find a version of the genome that's been chopped up into genes. It exists in the form of a .ffn file on the NCBI [ftp site](ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Fusobacterium_nucleatum_ATCC_25586_uid57885/). This file includes only the protein coding region of the genome segments. A quick analysis of this file shows that it's about .2Mbp shorter than the whole genome and is chopped into 1983 genes which are an average of 1kb long and the longest is 9kb which is below the 32 kb limit for velvet. 

Let's run the pipeline again:
~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/extract/FusoATCC

cp *.ffn ../Database
cd ../Database
mv NC_003454.ffn fuso.single.genes.fasta
bowtie2-build fuso.single.genes.fasta fuso.single.genes
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/velvet

bowtie2 /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single.genes -q /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/run2/DN/All.D1.Tongue.run2.norm.fq -p 16 -S velvet.fusogenes.sam  


cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/velvet

velveth k21_genes 21 -reference /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single.genes.fasta -short -sam velvet.fusogenes.sam

~~~~

Logfile:

~~~~

[0.000000] Reading FastA file /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database
/fuso.single.genes.fasta;
[0.066652] Overlapping reference coordinates:
[0.066654] gi|19703352|ref|NC_003454.1|:0-1104
[0.066655] gi|19703352|ref|NC_003454.1|:0-77
[0.066656] Exiting...
~~~~

So the ffn file isn't that useful because some of the genes are overlapping!!

