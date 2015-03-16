---
layout: post
title:  "March madness part 1"
date:   2015-3-8
comments: true
---

Pat challenged the lab to read as many papers as we can in the month of March. Here are the summaries that I also posted to the github page. 

##Wickham, H. Tidy data. *J. Stat. Software.* 2014 August; 59(10):1-23.
http://www.jstatsoft.org/v59/i10/paper

Overall: 7

Importance: 9

Innovation: 6

Who should read: everyone

Tags: R, statistics, data cleaning

[Hadley Wickham](http://en.wikipedia.org/wiki/Hadley_Wickham) is the current head of RStudio where he has been involved in developing packages for R that aid in data visualization (like ggplot2). In this paper he describes the importance, creation, and usage or *tidy* datasets. Tidy data is a way of formating a dataset in a table with each column representing a variable and each independent observation is a row. Wickham lists a few common ways that data sets can be tidied up, manipulated, and visualized. He doesn't go into too much detail of the *hows*, but he stresses the *whys*. Tidy data makes it easier to do efficient and repeatable analyses. Before doing an analysis, it's important to make sure your data set is tidy and know whether the tools you are using are input-tidy and output-tidy. Finally, he makes a call to improve tidy data tools.

I was instantly drawn to this paper, because I'm currently frustrated with all the translating of data that I have to do in my project. By creating a pipeline that uses lots of different software tools, I find myself writing simple scripts that parse, reformat, sort, extract data to funnel into the next part of the pipeline and then sometimes manually searching for typos in variable and sample names. So much time that I have wasted. I wish this paper had more diversity in the tool suggestions, but I did get a couple good ideas to look into. 

A couple R tools/packages that might be worth looking into:

1. Melting is used to expand a table into tidy data when variables are listed in rows and columns (these are often used in presentations to describe data). Melting can also be used to expand out data sets where a column contains two variables (such as sex_age). I was thinking this could be used to create and tidy up metadata from sample names. 
2. plyr package has tools to combine multiple tables with similar types of data (like the join() function which may be better than the merge() function that I already use). This is **so useful** for combining data from different experiments and avoiding matching up columns in excel like I did in college. 
3. ddply() which is an alternative to the by() function, but it is output-tidy (by() outputs a list).

##Time series community genomics analysis reveals rapid shifts in bacterial species, strains, and phage during infant gut colonization Sharon I, Morowitz MJ, Thomas BC, Costello EK, Relman DA, Banfield JF. Genome Res. 2013 Jan;23(1):111-20. doi: 10.1101/gr.142315.112. Epub 2012 Aug 30.

Overall: 9

Importance: 7

Innovation: 9

Who should read: @kdiverson, 

Tags: genomics, metagenomic assembly, microbiome development 

This paper from the Banfield lab is a great example of constructing complete genomes from metagenomic data. The goal was to extract complete genomes from time series samples from a relatively simple community - the infant gut. They developed a novel pipeline which took into account coverage and composition data. First, they assembled the metagenomic reads into contigs and extended the assembly by iteratively assembling based on coverage. This method makes the assumption that strains at a similar abundance level will have the same read coverage throughout the genome. Second, they separated the contigs into bins based on tetranucleotide composition using ESOM. From this binning step they generated 8 almost completely genomes and 9 draft genomes as well as 3 viral genomes. 

This article demonstrates a few benefits of assembling full genomes compared to using metagenomic data for community gene profiling or only sequencing 16S data. By constructing full length genomes from time-series samples, changes in abundance of specific strains could be monitored over the time course. This is a benefit to using metagenomic data instead of 16S, because here there is higher strain resolution, but 16S might not show strain differences. Also, in this experiment they were able to assemble viral genomes, which cannot be assessed with 16S sequencing. Finally, they were able to do comparative genomics between two strains of Propionibacterium. By comparing full length genomes, they could identify genes associated with metabolic advantages and persistence in the gut. 

Sharon et al. demonstrate an innovative method to extract full length genomes from metagenomic data. The supplementary methods are fairly detailed, but unfortunately not easily reproduced because they did not publish scripts or software. Their method is pretty similar to CONCOCT which also take into account coverage and composition data, and it might be interesting to compare the methods using samples with more complexity than the infant gut. 


##The interplay between microbiome dynamics and pathogen dynamics in a murine model of Clostridium difficile Infection.
Reeves AE, Theriot CM, Bergin IL, Huffnagle GB, Schloss PD, Young VB.
Gut Microbes. 2011 May-Jun;2(3):145-58. Epub 2011 May 1.
http://www.ncbi.nlm.nih.gov/pubmed/21804357

Overall: 8

Importance: 8

Innovation: 7

Who should read: C.diff people

Tags: C. difficile, colonization resistance

This paper from our microbiome group at Michigan does a good job demonstrating that colonization resistance to C. difficile is dependent on the ecology of the gut microbial communities. The authors treated mice with different antibiotic coctails and challenged with different doses of C. difficile. Depending on the antibiotic used and the dose of C. difficile, the mice had differences in infection and severity of disease. Mice that were treated with antibiotics and recovered from CDI had a microbiome similar to untreated animals, while animals that had severe colitis developed a more dysbiotic community enriched in Proteobacteria. Additionally, there was a dose dependent effect on mortality in mice treated with cefoperazone, suggesting that the indigenous community can recover and outcompete C. difficile at a low abundance. Therefore, the authors propose a balance between the ability of the microbiome to recover from antibiotics and the rate of C. difficile expansion in the dysbiotic community. 

I like that this paper is very clear in demonstrating a mechanism (colonization resistance), instead of getting bogged down in details of different bugs or communities. The experiments are hypothesis driven with carefully controlled variables. It is  interesting that some mice get sick and some clear the infection when treated with the same antibiotics even though they should have started with similar communities. The authors don't dwell too much on this effect, just the result of it. I noticed that some of the mice that didn't get sick (CDI well), still had a fair amount of cellular infiltrate in the colon. I wonder if that was just leftover from infection, or if it is a result of the antibiotic treatment. 













