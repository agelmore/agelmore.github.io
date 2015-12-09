---
layout: post
title:  "Creating a pangenome - updated"
date:   2015-11-29
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

##Format input

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

The output from the program is a list of the sequence IDs that clustered (`faa_homologues/NC003454_f0_alltaxa_algBDBH_e0_.cluster_list`) and a folder with all of the consensus sequences from those clusters (`faa_homologues/NC003454_f0_alltaxa_algBDBH_e0_`). This is great but it isn't really what I wanted. I want a file that has all the shared genes *AND* all the novel genes from each of the species. That way I have a list of ALL the Fusobacterium genes. 

##OMCL algorithm	

I ran the script again with the option to use the OrthoMCL algorithm. It used slightly more RAM, but the results are almost identical. It split the genes into 1474 clusters instead of 1469. 

##Adding novel genes

The number of taxa that have a sequence is set to default to the number of total genomes. This makes it so that every genome must have the homologues sequence for it to count as a cluster. This creates a *core* genome. What I am interested in is a *pan* genome which contains the core genes *and* the novel genes. To do this you have to set the number of taxa parameter (-t) to 0. 

~~~~
get_homologues.pl -d faa -t 0 -M 
~~~~

With these options it created the clusters of homologues genes (paralogs and orthologs) and then added the singletons that didn't have any match.

~~~~
# checking input files...
# NC_003454.faa 1983
# NC_021281.faa 2107
# NC_022196.faa 2151

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

I got distracted trying to figure out what these extra scripts do:

##compare_cluster.pl

The GET_HOMOLOGUES script *compare_cluster.pl*  will intersect the cluster files and make some graphs. It will create parsimony trees to look at divergence of different genes.

~~~~
compare_clusters.pl -o faa_compare -d faa_homologues/NC003454_f0_0taxa_algOMCL_e0_ -t 0 -m -T
~~~~

Summary of options:

-o	: directory for the output
-t 	: number of taxa that have to have the sequence to create a cluster
-m	: produces intersection pangenomic matrix
-T	: produces parsimony based pangenomic tree

I don't think this script is very useful to me. It's usually used when combining genomes from different samples to see if there is gene overlap. I could have used this for comparative genomics if I were able to assemble full genomes from my metagenomic samples. 

##parse pangenome

This script is supposed to parse the genes into core and shell genes (genes present in all genomes vs. genes present in one or two. Might be cool to visualize which genes are shared between the genomes I'm using and which are rare. 

~~~~
parse_pangenome_matrix.pl -m faa_homologues/NC003454_f0_0taxa_algOMCL_e0_Avg_identity.tab -s
~~~~

Can't get it to work. I could do this summary pretty easily myself though by using information from the `NC003454_f0_0taxa_algOMCL_e0_.cluster_list` file which tells me how many genomes contain genes from each cluster.

#Create pangenome

I would like to generate some sort of consensus sequence for clusters that have multiple genes, but for now I'm just going to pick one of them and use that as a representative sequence for the cluster. They all should have similar sequences, so hopefully any orthologs will align to the representative later when I'm aligning reads. 

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/faa_homologues/NC003454_f0_0taxa_algOMCL_e0_

for f in *.faa; do head -n2 $f | cat >> pan.complete.first.faa; done
~~~~


#Add draft genomes

Ncbi has a ton of Fusobacterium draft genomes. Get_homologues makes it easy to add more genomes to the input folder, so I'll just put them in there. I have to cat the draft scaffolds to a single file for each sample.

I put the address for all the drafts in a file called wget/draft_fuso.txt

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/draft

wget -i ../wget/draft_Fuso.txt

#unzip all the files and put in a folder
for f in *.tgz; do name=${f:0:7}; mkdir $name; quicksubmit "tar zxf $f -C $name" $quickpara; done

#move the zipped files to the same place
mkdir zipped
mv *.tgz zipped

#combine draft files into one per species
for i in NZ*; do cd $i; quicksubmit "cat *.faa >> ../../faa/$i.faa"; cd ..; done

#rerun get_homologues
get_homologues.pl -d faa -t 0 -M
~~~~

This finished in about 12 hours without using any of the parallel options. That's actually pretty quick because with 30 draft genomes and 3 complete genomes, the program blasts all against all for 1089 total blast files. 

Here are the genome summaries:

~~~~
# NC_003454.faa 1983
# NC_021281.faa 2107
# NC_022196.faa 2151
# NZ_AABF.faa 2250
# NZ_AARG.faa 2362
# NZ_ACDB.faa 2538
# NZ_ACDC.faa 2354
# NZ_ACDD.faa 1890
# NZ_ACDE.faa 2154
# NZ_ACDF.faa 2418
# NZ_ACDG.faa 1833
# NZ_ACDH.faa 3191
# NZ_ACDS.faa 2334
# NZ_ACET.faa 1621
# NZ_ACIE.faa 3008
# NZ_ACIF.faa 2624
# NZ_ACJY.faa 2568
# NZ_ACPU.faa 2131
# NZ_ACQE.faa 2123
# NZ_ACUO.faa 2641
# NZ_ADDB.faa 2423
# NZ_ADEE.faa 2037
# NZ_ADGF.faa 2006
# NZ_ADLZ.faa 3197
# NZ_ADVK.faa 2126
# NZ_AFQD.faa 2777
# NZ_AGAD.faa 2274
# NZ_AGEH.faa 2307
# NZ_AGWJ.faa 3755
# NZ_AJSY.faa 2060
# NZ_AKXI.faa 1821
# NZ_ALKK.faa 2173
# NZ_ALVD.faa 2200

# 33 genomes, 77437 sequences
~~~~

And the outfile for the number of clusters and runtime:

~~~~
# number_of_clusters = 15405
# cluster_list = faa_homologues/NZACET_f0_0taxa_algOMCL_e0_.cluster_list
# cluster_directory = faa_homologues/NZACET_f0_0taxa_algOMCL_e0_

# runtime: 36979 wallclock secs (1898.38 usr 16.82 sys + 41202.69 cusr 327.50 csys = 43445.39 CPU)
# RAM use: 1332.7 MB
~~~~

Lets look at some stats of the clusters. (Apparently you can do this with one of the get_homologues scripts but I can't get them to work.)

~~~~
for f in {1..100}; do grep -c "size=$f " NZACET_f0_0taxa_algOMCL_e0_.cluster_list >> genespercluster.txt; done

#in R
x<- read.delim(file="genespercluster.txt", header=F)
x$size=rownames(x)
plot(x$size, x$V1, type="h", log="y", main="Genes per cluster", xlab="number of gene copies", ylab='number of clusters')
~~~~


![Coverage map reads extracted from full database and assembled]({{ site.url }}/images/Genespercluster.png)

Here's a plot showing how the genes clustered together. **Notice the Y axis is on a log scale so I could see the separation better.**  Most of the clusters have fewer than 33 genes (one per taxa), but some have many more which means they are paralogs with multiple copies in the same genome. I checked in to see what genes are in the clusters with high gene copy number. 

~~~~
awk '$3== "size=94" {print}' NZACET_f0_0taxa_algOMCL_e0_.cluster_list

#output
cluster 58_NP_602381.1 size=94 taxa=24 file: 58_NP_602381.1.faa dnafile: void
~~~~

According to the annotations, it looks like this cluster is a bunch of outer membrane proteins and/or autotransporters. The cluster that has 64 copies are all "translation elongation factor G" and with 62 are mostly "tRNA modification enzyme". That makes sense that those would be paralogs. Cool!

**There are a lot of clusters that only have a single gene**, suggesting novel genes, high divergence (so that they didn't blast to eachother), or incomplete genes. I wonder if these are an artifact of the get_homologues method. Because I used so many draft genomes, it's possible that there are a bunch of incomplete genes which come up as unique, but aren't really. **How am I going to check that?** There are 15405 clusters total and 8889 of those are singletons. If I remove the singletons I have 6519 clusters, which is pretty similar to the number of clusters that they found in the GET_HOMOLOGUES papers when using 50 Strep genomes. I think it might be best to just remove the singletons. I'll miss any genes that are completely novel, but it will get rid of the incomplete gene error.

I'm running the program again with the option -t 1 which will remove singletons:

~~~~
get_homologues.pl -d faa -t 1 -M
~~~~

It finished. It worked exactly as I expected. The number of clusters is now 6516 which is in the range I would want. 

~~~~
# number_of_clusters = 6516
# cluster_list = faa_homologues/NZACET_f0_1taxa_algOMCL_e0_.cluster_list
# cluster_directory = faa_homologues/NZACET_f0_1taxa_algOMCL_e0_

# runtime: 42301 wallclock secs (1997.89 usr 16.16 sys + 44823.06 cusr 1782.71 c
sys = 48619.82 CPU)
# RAM use: 1332.4 MB
~~~~


#Create nucleotide pangenome

Run again with the .ffn files from ncbi. These are the same as the .faa but as DNA sequences instead of protein. For get_homologues you have to include the protein sequences with the fasta sequences in the same order for it to give you a nucleotide output (I guess it will cluster the protein sequences and then grab the nucleotide sequences later). I discovered through trial and error that the .ffn files have to have the .fna extension even though they are .ffn. And they have to have the same name as the .faa files. If I set this all up perfectly, nucleotide clusters will be produced. 

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/ffn

wget -i ../wget/complete_all_ffn.txt
wget -i ../wget/draft_ffn.txt

#unzip all the files and put in a folder
for f in *ffn.tgz; do name=${f:0:7}; mkdir $name; quicksubmit "tar zxf $f -C $name" $quickpara; done

#move the zipped files to the same place
mkdir zipped
mv *.tgz zipped

#combine draft files into one per species
for i in NZ*; do cd $i; quicksubmit "cat *.ffn >> ../../$i.fna"; cd ..; done

#get the matching faa files that I used before
cp ../faa/NZ*.faa .

#rerun get_homologues
get_homologues.pl -d ffn -t 1 -M
~~~~

It worked. In the cluster list file are clusters in both protein and nucleotide format (matching names, only difference is ".faa" vs ".fna". The genes are the same, just presented in each format. This is good because now I can run BWA (which only works on nucleotide databases). I moved the nucleotide clusters to their own folder `/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/ffn_homologues/NZACET_f0_1taxa_algOMCL_e0_/nucleotide
`

#Do it again with different t

**Started Nov 29, 2015 at 3pm**


Run get_homologues with different t cutoff. This way I will be able to choose which pangenome is the best. As I increase t I through away more rare (or incomplete genes), but I also get a cleaner dataset to work with. 

~~~~
for f in {1,2,3,4,5,6,7,8,9,10,15,20,25,30,33}; do cp -r ffn ffn_t_$f; quicksubmit "get_homologues.pl -d ffn_t_$f -t $f -M; rm -r ffn_t_$f" $quickpara; done
~~~~

I will make a rarefaction curve and SCG analysis for each clustering cutoff to help me decide which is the best to use. 

##Rarefaction curves

~~~~
for f in {1,2,3,4,5,6,7,8,9,10,15,20,25,30,33}; do quicksubmit "compare_clusters.pl -o rarefaction_t$f -d ffn_t_"$f"_homologues/*_ -t 0 -m; rm rarefaction_t"$f"/*faa" $quickpara; done

~~~~

I need to figure out how to turn the output tab file into a shared file for mothur to use.


##Make index files for vart

##Combine genes into single file
~~~~
for i in {1,2,3,4,5,6,7,8,9,10,15,20,25,30,33}; do f=/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/vart; cd $f/ffn_t_"$i"_homologues; cd *_; quicksubmit "cat *.fna >> $f/index/all.t"$i".fna" $quickpara; done
~~~~

Number of genes per file. As t increases, the number of genes that were used to make the pangenome goes down, because if they weren't present in the set number of genomes, those genes were discarded.

~~~~
all.t10.fna:52867
all.t15.fna:47987
all.t1.fna:68548
all.t20.fna:44023
all.t25.fna:37166
all.t2.fna:1091648
all.t30.fna:29042
all.t33.fna:9869
all.t3.fna:64897
all.t4.fna:61409
all.t5.fna:58674
all.t6.fna:57243
all.t7.fna:55785
all.t8.fna:54735
all.t9.fna:53739
~~~~

##Create index file for each t

~~~~
#create a "group" file. Every line is a unique sequence and the cluster that it belongs to. I need this for when I'm doing my alignment.

dir=/mnt/EXT/Schloss-data/amanda/Fuso/pangenome/vart

Doesn't work quicksubmit for some reason:

for f in *.fna; do cat $f | awk 'NR % 2 ==1' | awk -v var=$f '{gsub(/>/,"")} {OFS="\t"} {print $1, var}' >> $dir/index/t"$a".index; done

~~~~



