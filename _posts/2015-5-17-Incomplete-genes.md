---
layout: post
title:  "Incomplete genes in draft genomes"
date:   2015-5-17
comments: true
---

When I set get_homologues to allow clusters with a single gene, I get a ton of singleton clusters (58% of the total number of clusters). I think that's because a lot of the genes in the draft genomes are incomplete. For now I removed those clusters, but I want to check to see if those incomplete genes are in the other clusters.

**Task:** Blast the singleton clusters against the other clusters to see if they are present there (and incomplete) or if they are in fact novel genes.

~~~~
#run get_homologues again with -t 0
mv ffn ffn_t0
get_homologues.pl -d ffn_t0 -t 0 -M

#merge clusters t>0 into single file. 
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/ffn_homologues/NZACET_f0_1taxa_algOMCL_e0_/nucleotide
cat *.fna >> genes_t1.fna

#merge clusters t=0 into single file
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/ffn_t0_homologues/NZACET_f0_0taxa_algOMCL_e0_
grep -c '>' *.fna > genes_per_cluster.txt
awk -F : '$2==1 {print $1}' genes_per_cluster.txt > genes_t0_id.txt
for f in $(cat genes_t0_id.txt); do cat $f >> genes_t0.fna; done
cp genes_t0.fna /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/incomplete

#index database
makeblastdb -in genes_t1.fna -dbtype nucl -out genes_t1.blast

#run blast
cd /mnt/EXT/Schloss-data/amanda/Fuso/pangenome/incomplete
blastn -db ../ffn_homologues/NZACET_f0_1taxa_algOMCL_e0_/nucleotide/genes_t1.blast -query genes_t0.fna -out incomplete.blast.1 -evalue 1e-5 -outfmt 6 -num_threads 16 -max_target_seqs 1
~~~~

**2646** of the singletons hit a gene in the clusters with an evalue of less than 1e-5. There were **8889** total singletons. I'm going to try bwa with all the singletons included and see if I get a lot of hits to those genes.