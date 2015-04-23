---
layout: post
title:  "SCG evaluation, one cluster?"
date:   2015-4-12
comments: true
---

The megahit assembly of reads that aligned with Fusobacterium genomes finished, and I'm interested to see how pure the assembly is. Because I already figured out how to use CONCOCT's [single copy core gene validation](http://agelmore.github.io/2015/02/28/Single-copycore.html), I wondered if I could use this pipeline **without clustering**. 

Run prodigal. I'm calling the contigs ref.bowtie.megahit. Don't forget the silly bit of editing I have to do to the contig ids. Also, the script requires the CONCOCT output with cluster numbers. Since I haven't ran CONCOCT, I just made a csv file that has all the contig names with the cluster #1. Clever, huh?

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/megahit/SCG
cp ../final.contigs.fa ref.bowtie.megahit.fa

#change symbols to help things later
sed -i 's/_/-/g' ref.bowtie.megahit.fa 

#create file with contig ids and cluster#1. 
grep '>' ref.bowtie.megahit.fa | sed 's/>//g' | sed 's/$/,1/' > cluster1.csv 


#now run prodigal
cd /mnt/EXT/Schloss-data/amanda/prodigal/prodigal.v2_50

./prodigal -a /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/megahit/SCG/ref.bowtie.megahit.faa -i /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/megahit/SCG/ref.bowtie.megahit.fa -f gff -p meta -o /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/megahit/SCG/ref.bowtie.megahit.gff 
~~~~

Prodigal creates two files. The .gff file has gene annotations for genes on the contigs input. The .faa file translates these genes and creates a file with each of the genes and their amino acid translations. Now run RPSBLAST.sh script from CONCOCT.


~~~~
concoctenv
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/megahit/SCG
$CONCOCT/scripts/RPSBLAST.sh -f ref.bowtie.megahit.faa -p -c 8 -r 1 

#have to do more editing so the script will work
sed -e 's/_/-/2' ref.out > ref.cogtable.out
~~~~

Now, generate the table using CONCOCT script COG_table.py

~~~~
python $CONCOCT/scripts/COG_table.py -b ref.cogtable.out -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c cluster1.csv --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > ref.bowtie.megahit_scg.tsv
~~~~


Graph the output:

~~~~
Rscript $CONCOCT/scripts/COGPlot.R -s ref.bowtie.megahit_scg.tsv -o ref.bowtie.megahit_scg.pdf
~~~~

![Single-copy core genes]({{ site.url }}/images/ref.bowtie.megahit_scg.png)

I got the script to work. Doesn't look like we have a single complete genome here. The two that have a lot of repeats are ribosomal genes which could be repeats in the same strain or different strain. Unfortunately, some of the missing cogs should be in Fusobacterium. I'm going to do the pipeline on a completely sequenced genome to see how many COGs are present in a genome we know is complete and Fuso.

##SCG control on completed Fuso genome

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/megahit/SCG
cp /mnt/EXT/Schloss-data/amanda/Fuso/extract/Database/fuso.single.fna fuso.control.fa

#change symbols to help things later
sed -i 's/_/-/g' fuso.control.fa 

#create file with contig ids and cluster#1. It's just one sequence
nano fuso.control.csv 
gi|19703352|ref|NC-003454.1|,1

#now run prodigal
cd /mnt/EXT/Schloss-data/amanda/prodigal/prodigal.v2_50

./prodigal -a /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/megahit/SCG/fuso.control.faa -i /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/megahit/SCG/fuso.control.fa -f gff -p meta -o /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/megahit/SCG/fuso.control.gff 

concoctenv
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/reference/bowtie/megahit/SCG
$CONCOCT/scripts/RPSBLAST.sh -f fuso.control.faa -p -c 8 -r 1 

#have to do more editing so the script will work
sed -e 's/_/-/2' fuso.out > fuso.control.cogtable.out

python $CONCOCT/scripts/COG_table.py -b fuso.control.cogtable.out -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c fuso.control.csv --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > fuso.control_scg.tsv

Rscript $CONCOCT/scripts/COGPlot.R -s fuso.control_scg.tsv -o fuso.control_scg.pdf
~~~~


![Single-copy core genes control]({{ site.url }}/images/fuso.control_scg.png)


This pipeline seems to be a pretty good control to see how complete my assembly is! The complete Fuso genome from NCBI had almost all of the COGs except for COG0198 and COG0504.
