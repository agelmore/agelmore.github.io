---
layout: post
title:  "Validate using Single-copy core genes"
date:   2015-2-28
comments: true
---

The CONCOCT pipeline has a script to evaluate the purity of clusters by using single-copy core genes. The idea is that if each cluster is a pure genome, then it will contain one and only one copy of core genes (i.e. 16s). 

The pipeline has you annotate and translate your contigs using prodigal, use the script RPSBLAST.sh to cog annotate the genes, 

Downloaded [prodigal](http://prodigal.ornl.gov/downloads.php) and ran the gene annotations on the megahit assembly of contigs greater than 1000bp. 

~~~~
./prodigal -a $HMP/D1.tongue/run2/concoct/1kb/annotations/proteins/megahit.1000.contigs_c10K.faa -i $HMP/D1.tongue/run2/megahit/megahit.1000.contigs_c10K.fa -f gff -p meta -o $HMP/D1.tongue/run2/concoct/1kb/annotations/proteins/megahit.1000.contigs_c10K.gff 
~~~~

Prodigal creates two files. The .gff file has gene annotations for genes on the contigs input. The .faa file translates these genes and creates a file with each of the genes and their amino acid translations. 

These protein sequences are now run against the cog database using the CONCOCT script RBSBLAST.sh. Downloaded Cog database (had to get one that was already formatted for rpsblast) from ncbi (wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz) and extracted. **I wish they had put this information in the read the docs, would have made it a lot easier to find. I spent a lot of time trying to figure out how to format the COG database for rpsblast. formatrpsdb no longer exists in the BLAST package.**

~~~~
$CONCOCT/scripts/RPSBLAST.sh -f megahit.1000.contigs_c10K.faa -p -c 8 -r 1
~~~~

Alternative:

~~~~
rpsblast -query $HMP/D1.tongue/run2/concoct/1kb/annotations/proteins/megahit.1000.contigs_c10K.faa -db Cog -evalue 0.00001 -outfmt "6 qseqid sseqid evalue pident score qstart qend sstart send length slen" -out blast_output.out
~~~~



This script took a little while. The output was put in a file called megahit.out. I moved this file to a folder called cog-annotations.

Next we use the CONCOCT script COG_table.py. This script requires the file `scg_cogs_min0.97_max1.03_unique_genera.txt` which contains the single-copy cog ids to use and the file `$CONCOCT/scgs/cdd_to_cog.tsv` which matches cog ids to their conserved domain database number. 

~~~~
python $CONCOCT/scripts/COG_table.py -b annotations/cog-annotations/megahit.out -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c concoct-output/clustering_gt1000.csv --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > evaluation-output/2clustering_gt1000_scg.tsv
~~~~

It worked and created the cog table, but it says that none of the contigs matched the COGs at all! I realized the problem is in the name of my contigs. Prodigal adds a _ character to the end of each contig when it cut the contigs into genes. The COG_table.py program knows to separate the contig names from the gene name using the _ symbol, but my contig names already contain that separator. I need to go back to before I run prodigal and rename all the contigs. I will also need to rename them in the cluster file. Ugg.  



~~~~
sed -i 's/_/-/g' megahit.1000.contigs_c10K.fa 
sed -i 's/_/-/g' clustering_gt1000.csv 

./prodigal -a $HMP/D1.tongue/run2/concoct/1kb/cogtable/megahit.1000.contigs_c10K.faa -i $HMP/D1.tongue/run2/concoct/1kb/cogtable/megahit.1000.contigs_c10K.fa -f gff -p meta -o $HMP/D1.tongue/run2/concoct/1kb/cogtable/megahit.1000.contigs_c10K.gff 

$CONCOCT/scripts/RPSBLAST.sh -f megahit.1000.contigs_c10K.faa -p -c 8 -r 1

sed -e 's/_/-/2' megahit.out > megahit.cogtable.out #replaces every second instance of _ in each line so that the COG_table.py script will work

python $CONCOCT/scripts/COG_table.py -b megahit.cogtable.out -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c clustering_gt1000.csv --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > clustering_gt1000_scg.tsv

~~~~

It worked!!!

Graph the output:

~~~~
Rscript $CONCOCT/scripts/COGPlot.R -s clustering_gt1000_scg.tsv -o clustering_gt1000_scg.pdf
~~~~


![Single-copy core genes heat map]({{ site.url }}/images/clustering_gt1000_scg.png)

