---
layout: post
title:  "Creating a pangenome"
date:   2015-5-3
comments: true
---


A pangenome of a species is the sum of the core genes found in all strains and strain-specific genes found in some but not all strains [(review)](http://www-ncbi-nlm-nih-gov.proxy.lib.umich.edu/pubmed/19086349).

I will be using the software package [GET_HOMOLOGUES](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3837814/) (Bruno Contreras-Moreira and Vinuesa, 2013. *Genome Research*). 

![GET_HOMOLOGUES Flowchart ]({{ site.url }}/images/get_homologues_flowchart.png)

The software package is open-source, automatic, customizable, and can handle large jobs. All the scripts are written in perl, yay. Basically, genes input in FASTA format are BLASTed against each other. Based on BLAST hits, genes are clustered into orthologous groups using one of three popular clustering algorithms: COGtriangles, Bidirectional Best hit, and OrthoMCL (recommended by Evan Snitkin). The final result is a fasta file with all the pangenomic genes (other formats and statistics optional). The paper also showed that the pipeline works for draft genomes, which is good because I have a lot of them.

##Install

~~~~
cd /mnt/EXT/Schloss-data/amanda/get_homologues

wget http://maya.ccg.unam.mx/soft/get_homologues-x86_64-20150306.tgz
tar xvfz get_homologues-x86_64-20150306.tgz
./install.pl

#test it
./get_homologues.pl -d sample_buch_fasta 
~~~~

##Format FASTA

The best way to input files into GET_HOMOLOGUES is by putting all of the sequences in a single folder in FASTA format. Each genome is in a separate file. This way, if more genomes are added later (or if we sequence some in lab), they can be easily added to the pangenome. 

I'm going to start with a pangenome of just F. nucleatum using complete genomes. See if that works and I can add the draft genomes later. I'm downloading all the ncbi entries into a single folder.

I tried to run the program using the fasta files from ncbi, but those haven't been annotated like the faa (amino acid) files have. I'm going to try running get_homologues with the amino acid files and the GenBank files. This will produce a protein sequences pangenome, but maybe I can translate it back later? 

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_gbk

wget -i complete_nucleatum_faa.txt
wget -i complete_nucleatum_gbk.txt

~~~~

##Run

~~~~
get_homologues.pl -i faa_gbk
~~~~

It didn't put anything into the outfile. Here's the log:

~~~~
# /mnt/EXT/Schloss-data/amanda/get_homologues/get_homologues-x86_64-20150306/get
_homologues.pl -i 0 -d faa_gbk -o 0 -e 0 -f 0 -r 0 -t all -c 0 -I 0 -m local -n 
2 -M 0 -G 0 -P 0 -C 75 -S 1 -E 1e-05 -F 1.5 -N 0 -B 50 -b 0 -s 0 -D 0 -g 0 -a '0
' -x 0 -R 0 -A 0

# results_directory=/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_gbk_homologu
es
# parameters: MAXEVALUEBLASTSEARCH=0.01 MAXPFAMSEQS=250 BATCHSIZE=100 KEEPSCNDHS
PS=0

# checking input files...
# NC_003454.faa 1983
# NC_003454.gbk 1983
# NC_009506.faa 11
# NC_009506.gbk 11
# NC_022196.faa 2151
# NC_022196.gbk 2151

# 6 genomes, 8290 sequences

# taxa considered = 6 sequences = 8290 residues = 2643340 MIN_BITSCORE_SIM = 18.
1

# mask=NC009506_f0_alltaxa_algBDBH_e0_ (_algBDBH)

# running makeblastdb with /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_gbk_h
omologues/NC_003454.faa.fasta

# running makeblastdb with /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_gbk_h
omologues/NC_003454.gbk.fasta

# running makeblastdb with /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_gbk_h
omologues/NC_009506.faa.fasta

# running makeblastdb with /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_gbk_h
omologues/NC_009506.gbk.fasta

# running makeblastdb with /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_gbk_h
omologues/NC_022196.faa.fasta

# running makeblastdb with /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_gbk_h
omologues/NC_022196.gbk.fasta

# running BLAST searches ...
# done

# concatenating and sorting blast results...
# sorting _NC_003454.faa results (1.2MB)
# sorting _NC_003454.gbk results (1.3MB)
# sorting _NC_009506.faa results (0.0017MB)
# sorting _NC_009506.gbk results (0.0017MB)
# sorting _NC_022196.faa results (1.3MB)
# sorting _NC_022196.gbk results (1.3MB)
# done


# parsing blast result! (/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_gbk_hom
ologues/tmp/all.blast , 5.1MB)
# parsing file finished

# creating indexes, this might take some time (lines=1.09e+05) ...

# construct_taxa_indexes: number of taxa found = 6
# number of file addresses = 1.1e+05 number of BLAST queries  = 8.3e+03

# clustering orthologous sequences

# clustering inparalogues in NC_009506.faa (reference)
# 0 sequences

# clustering inparalogues in NC_003454.faa
# 54 sequences

# finding BDBHs between NC_009506.faa and NC_003454.faa
# 0 sequences

# clustering inparalogues in NC_003454.gbk
# 54 sequences

# finding BDBHs between NC_009506.faa and NC_003454.gbk
# 0 sequences

# clustering inparalogues in NC_009506.gbk
# 0 sequences

# finding BDBHs between NC_009506.faa and NC_009506.gbk
# 11 sequences

# clustering inparalogues in NC_022196.faa
# 30 sequences

# finding BDBHs between NC_009506.faa and NC_022196.faa
# 1 sequences

# clustering inparalogues in NC_022196.gbk
# 30 sequences

# finding BDBHs between NC_009506.faa and NC_022196.gbk
# 1 sequences

# looking for valid ORF clusters (n_of_taxa=6)...


# number_of_clusters = 0

# runtime: 363 wallclock secs (17.07 usr  0.32 sys + 634.36 cusr  5.79 csys = 65
7.54 CPU)
# RAM use: 116.9 MB
qsub working directory absolute is
/mnt/EXT/Schloss-data/amanda/Fuso/pangenome
~~~~







