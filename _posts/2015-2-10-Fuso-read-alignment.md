---
layout: post
title:  "How many Fuso reads?"
date:   2015-2-10
comments: true
---

I made a quick pipeline to see how many Fuso reads were used in the megahit assembly. Basically, I made a list of contigs that mapped to Fuso using blast, used a script to pull out those contigs from the assembly file, and then ran bowtie2 to find the percent reads mapping. I did this with the complete list of Fuso and with the Fuso contigs > 1kb. 

*Output:*

5.17% total reads map to fuso contigs

1.12% total reads map to fuso contigs greater than 1kb

~~~~
##run concoct script on all contigs > 1kb from the assembly
concoctenv
cd $HMP/D1.tongue/run2
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m megahit/final.contigs.1000.fa > megahit/megahit.1000.contigs_c10K.fa

##blast contigs > 1kb against fuso database
blastn -db Fuso.all.db.make -query /mnt/EXT/Schloss-data/amanda/Fuso/HMP/D1.tongue/run2/megahit/megahit.1000.contigs_c10K.fa -out blast.fuso.1000 -evalue 1e-5 -outfmt 6 -num_threads 16 -max_target_seqs 1

##cut out only the contig names
cat blast.fuso.1000 | cut -f1 > blast.fuso.1000.contigsonly

Found a script online that can pull out sequences from a list of contigs

#change the output format of the fasta file
sed 's/>/\'$'\n>/g' fuso.all.fa > fuso.all.2.fa
sed 's/>/\'$'\n>/g' fuso.1000.fa > fuso.1000.2.fa 


bowtie2-build blast/fuso.all.fa blast/fuso.all.fa
bowtie2-build blast/fuso.1000.fa blast/fuso.1000.fa


bowtie2 blast/fuso.all.fa -q cat/All.D1.tongue.run2.cat.fq -p 16 -S readstofuso.all.sam 
bowtie2 blast/fuso.1000.fa -q cat/All.D1.tongue.run2.cat.fq -p 16 -S readstofuso.1000.sam 

~~~~

