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

It loaded the files just fine, but then it didn't find ANY matches between the genomes. Obviously, thats not possible since they are all F. nucleatum. Maybe the proteins aren't annotated correctly? Or maybe blast is being too stringent? 

Oh I just realized that NC_009506.fna is a plasmid and not a full genome. That was my problem. I reran with only NC_003454.faa  NC_021281.faa  NC_022196.faa from NCBI. NC_021281.faa is not a nucleatum genome, but I decided to try making a genus level genome because I have so few complete genomes. 

I reran the script with the three complete genomes in amino acid sequence format.

~~~~
get_homologues.pl -d faa
~~~~

The script took about 10 minutes. 

Outfile:

~~~~
# results_directory=/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_homologues
# parameters: MAXEVALUEBLASTSEARCH=0.01 MAXPFAMSEQS=250 BATCHSIZE=100 KEEPSCNDHSPS=0

# checking input files...
# NC_003454.faa 1983
# NC_021281.faa 2107
# NC_022196.faa 2151

# 3 genomes, 6241 sequences

# taxa considered = 3 sequences = 6241 residues = 1987852 MIN_BITSCORE_SIM = 17.9

# mask=NC003454_f0_alltaxa_algBDBH_e0_ (_algBDBH)

# running makeblastdb with /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_homologues/NC_003454.faa.fasta

# running makeblastdb with /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_homologues/NC_021281.faa.fasta

# running makeblastdb with /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_homologues/NC_022196.faa.fasta

# running BLAST searches ...
# done

# concatenating and sorting blast results...
# sorting _NC_003454.faa results (0.9MB)
# sorting _NC_021281.faa results (0.9MB)
# sorting _NC_022196.faa results (0.95MB)
# done


# parsing blast result! (/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_homologues/tmp/all.blast , 2.7MB)
# parsing file finished

# creating indexes, this might take some time (lines=5.88e+04) ...

# construct_taxa_indexes: number of taxa found = 3
# number of file addresses = 5.9e+04 number of BLAST queries  = 6.2e+03

# clustering orthologous sequences

# clustering inparalogues in NC_003454.faa (reference)
# 81 sequences

# clustering inparalogues in NC_021281.faa
# 94 sequences

# finding BDBHs between NC_003454.faa and NC_021281.faa
# 1533 sequences

# clustering inparalogues in NC_022196.faa
# 70 sequences

# finding BDBHs between NC_003454.faa and NC_022196.faa
# 1594 sequences

# looking for valid ORF clusters (n_of_taxa=3)...


# number_of_clusters = 1469
# cluster_list = faa_homologues/NC003454_f0_alltaxa_algBDBH_e0_.cluster_list
# cluster_directory = faa_homologues/NC003454_f0_alltaxa_algBDBH_e0_

# runtime: 309 wallclock secs ( 4.27 usr  0.28 sys + 559.54 cusr 11.20 csys = 575.29 CPU)
# RAM use: 37.4 MB
~~~~

The output from the program is a list of the sequence IDs that clustered (`faa_homologues/NC003454_f0_alltaxa_algBDBH_e0_.cluster_list`) and a folder with all of the consensus sequences from those clusters (`faa_homologues/NC003454_f0_alltaxa_algBDBH_e0_`). This is great but it isn't really what I wanted. I want a file that has all the shared genes *AND* all the novel genes from each of the species. That way I have a list of ALL the Fusobacterium genes. 

