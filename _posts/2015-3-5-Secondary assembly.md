---
layout: post
title:  "Secondary assembly"
date:   2015-3-9
comments: true
---

Now that I have identified two pure (ish) clusters that contain Fusobacterium, I am going to reassemble the reads used to make those contigs and (hopefully) get some long Fuso contigs out. 

#Cluster contig files

Pull out contigs into separate files for each cluster.

~~~~
cd /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/run2/concoct/1kb/assembly2
awk -F , '$2 == "17"' ../concoct-output/clustering_gt1000.csv > cluster.17 
cut -d "," -f1 cluster.17 > cluster.17.cut
awk -F , '$2 == "21"' ../concoct-output/clustering_gt1000.csv > cluster.21 
cut -d "," -f1 cluster.21 > cluster.21.cut
cat cluster.17.cut cluster.21.cut > cluster.17.21
~~~~

Use mothur command to extract sequences in the id file and create a single fasta with all the contigs from clusters 17 and 21. 

~~~~
mothur '#get.seqs(accnos=cluster.17.21, fasta=megahit.1000.contigs_c10K.fa)'
mv megahit.1000.contigs_c10K.pick.fa cluster.17.21.fa
~~~~


Run bowtie to generate a sam file of all the reads used to create these contigs

~~~~
bowtie2-build cluster.17.21.fa cluster.17.21
bowtie2 cluster.17.21 -q $HMP/D1.tongue/run2/cat/All.D1.tongue.run2.cat.fq -p 16 -S cluster.17.21.sam 
~~~~

Outfile statistics:

~~~~
1622964606 reads; of these:
  1622964606 (100.00%) were unpaired; of these:
    1612764003 (99.37%) aligned 0 times
    9689093 (0.60%) aligned exactly 1 time
    511510 (0.03%) aligned >1 times
0.63% overall alignment rate
~~~~

This took me a couple tries to figure out how to do this. There is probably a better way, but this worked. 

~~~~
#make BAM file
samtools view -bT cluster.17.21 cluster.17.21.sam > cluster.17.21.bam

#use the -F4 option. The -F option removes the specified FLAG. The 4 flag is unmapped reads. 
samtools view -F4 cluster.17.21.bam > cluster.17.21.mapped.sam.2

#cut out name, sequence, and quality from sam file
cut -f1,10,11 cluster.17.21.mapped.sam > cluster.17.21.mapped.cut.sam
~~~~

Make into fastq format using python script samtofastq.py 

	-is there a samtool or picard to do this? - SamToFastq threw an error that it couldn't parse the sam file after running the F4 option

~~~~
#!/usr/bin/python
from __future__ import print_function
import sys

tab=open('cluster.17.21.mapped.cut.sam','r') #fastq sequences separated by tab
fastq=open('cluster.17.21.mapped.cut.fastq','wt') #fastq format

seqid=[] #list of sequence ID
seq=[]   #list of read sequence
qual=[]  #list of quality score
	
for row in tab:
	row=row.strip().split('\t')
	seqid.append(row[0]) 
	seq.append(row[1])  
	qual.append(row[2])
tab.close()

for i in range(len(seqid)):
	print("@", seqid[i], '\n', seq[i], '\n', '+', '\n', qual[i], sep='', file=fastq)
fastq.close()
~~~~

#Assemble!:

This fastq file has 10200603 reads which is about 1% of the total reads from the 20 samples. 

*Should I do digital normalization first? The fastq file is only 2.3G, but it might be good to normalize?*

~~~~
python ./megahit -m 45e9 -r $HMP/D1.tongue/run2/concoct/1kb/assembly2/cluster.17.21.mapped.cut.fastq --cpu-only -l 101 -o $HMP/D1.tongue/run2/concoct/1kb/assembly2/megahit
~~~~

The first time it finished with an error that the reads are longer than 100bp. That's weird because I assembled the first time with the read length at 100...but now some are 101bp? Well, anyways it's assembling. Hopefully it finishes tomorrow. 






