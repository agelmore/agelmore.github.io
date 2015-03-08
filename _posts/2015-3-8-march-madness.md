---
layout: post
title:  "March madness part 1"
date:   2015-3-8
comments: true
---

Pat challenged the lab to read as many papers as we can in the month of March. Here are the summaries that I also posted to the github page. 

Wickham, H. Tidy data. *J. Stat. Software.* 2014 August; 59(10):1-23.

Overall: 7

Importance: 9

Innovation: 6

Who should read: everyone

Tags: R, statistics, 

[Hadley Wickham](http://en.wikipedia.org/wiki/Hadley_Wickham) is the current head of RStudio where he has been involved in developing packages for R that aid in data visualization (like ggplot2). In this paper he describes the importance, creation, and usage or *tidy* datasets. Tidy data is a way of formating a dataset in a table with each column representing a variable and each independent observation is a row. Wickham lists a few common ways that data sets can be tidied up, manipulated, and visualized. He doesn't go into too much detail of the *hows*, but he stresses the *whys*. Tidy data makes it easier to do efficient and repeatable analyses. Before doing an analysis, it's important to make sure your data set is tidy and know whether the tools you are using are input-tidy and output-tidy. Finally, he makes a call to improve tidy data tools.

I was instantly drawn to this paper, because I'm currently frustrated with all the translating of data that I have to do in my project. By creating a pipeline that uses lots of different software tools, I find myself writing simple scripts that parse, reformat, sort, extract data to funnel into the next part of the pipeline and then sometimes manually searching for typos in variable and sample names. So much time that I have wasted. I wish this paper had more diversity in the tool suggestions, but I did get a couple good ideas to look into. 

A couple R tools/packages that might be worth looking into:

1. Melting is used to expand a table into tidy data when variables are listed in rows and columns (these are often used in presentations to describe data). Melting can also be used to expand out data sets where a column contains two variables (such as sex_age). I was thinking this could be used to create and tidy up metadata from sample names. 
2. plyr package has tools to combine multiple tables with similar types of data (like the join() function which may be better than the merge() function that I already use). This is **so useful** for combining data from different experiments and avoiding matching up columns in excel like I did in college. 
3. ddply() which is an alternative to the by() function, but it is output-tidy (by() outputs a list).
