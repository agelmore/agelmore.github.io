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

##OMCL algorithm	

I ran the script again with the option to use the OrthoMCL algorithm. It used slightly more RAM, but the results are almost identical. It split the genes into 1474 clusters instead of 1469. 

##Adding novel genes

The number of taxa that have a sequence is set to default to the number of total genomes. This makes it so that every genome must have the homologues sequence for it to count as a cluster. This creates a *core* genome. What I am interested in is a *pan* genome which contains the core genes *and* the novel genes. To do this you have to set the number of taxa parameter (-t) to 0. 

~~~~
get_homologues.pl -d faa -t 0 -M -A
~~~~

With these options it created the clusters of homologues genes (paralogs and orthologs) and then added the singletons that didn't have any match.

~~~~

# identifying orthologs between NC_003454.faa and NC_021281.faa (0)
# 1579 sequences

# identifying orthologs between NC_003454.faa and NC_022196.faa (0)
# 1647 sequences

# identifying inparalogs in NC_003454.faa
# 81 sequences

# identifying orthologs between NC_021281.faa and NC_022196.faa (0)
# 1670 sequences

# identifying inparalogs in NC_021281.faa
# 94 sequences

# identifying inparalogs in NC_022196.faa
# 70 sequences

# running MCL (inflation=1.5) ...
# running MCL finished

# find_OMCL_clusters: parsing clusters (/mnt/EXT/Schloss-data/amanda/Fuso/pangen
ome/faa_homologues/tmp/all_ortho.mcl)

# add_unmatched_singletons : 1029 sequences, 3 taxa

# looking for valid ORF clusters (n_of_taxa=0)...

# number_of_clusters = 2838
# cluster_list = faa_homologues/NC003454_f0_0taxa_algOMCL_e0_.cluster_list
# cluster_directory = faa_homologues/NC003454_f0_0taxa_algOMCL_e0_

# runtime: 187 wallclock secs ( 3.36 usr  0.28 sys + 334.93 cusr  2.44 csys = 34
1.01 CPU)
# RAM use: 40.4 MB
~~~~

So now there are 2838 clusters that can be present in 1, 2, or 3 genomes as well as genes that are present multiple times in the same genome. The folder `faa_homologues/NC003454_f0_0taxa_algOMCL_e0_` contains the amino acid sequences for each of them separately, but what I want is a single file with a single sequence for each cluster. Also, I need them to be in fasta format.

##compare_cluster.pl

Now we use the GET_HOMOLOGUES script *compare_cluster.pl* which will intersect the cluster files and find a consensus sequence.

~~~~
compare_clusters.pl -o faa_compare -d faa_homologues/NC003454_f0_0taxa_algOMCL_e0_ -t 0 -m -T
~~~~

Summary of options:

-o	: directory for the output
-t 	: number of taxa that have to have the sequence to create a cluster
-m	: produces intersection pangenomic matrix
-T	: produces parsimony based pangenomic tree

Actually I don't think this script is very useful to me. It's usually used when combining genomes from different samples to see if there is gene overlap. I could have used this for comparative genomics if I were able to assemble full genomes from my metagenomic samples. 

##parse pangenome

This script is supposed to parse the genes into core and shell genes (genes present in all genomes vs. genes present in one or two. 

~~~~
parse_pangenome_matrix.pl -m faa_homologues/NC003454_f0_0taxa_algOMCL_e0_Avg_identity.tab -s
~~~~