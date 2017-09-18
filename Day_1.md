
Microbiome course: Day 1
=========================

This morning you have heard all the advantages of microbial study using the new technologies based on sequencing. The objective of the next two practice sessions is to learn some of the methods we use to investigate the characteristics of microbial communities. To do so, we will use faecal samples from a [group of volunteers][1] (controls) and a group of [patients with Crohn's disease][2] (cases) from the UMCG. 


Module 1: From sequencing reads to taxonomies
---------------------------------------------

Fortunately for you, during this course you will not have to deal or isolate the DNA from any faecal samples. After sequencing, we usually remove low quality reads and remove those sequences that belong / align to the host genome (in our case: the human genome). Once we have cleaned our sequenced reads, we can proceed to identify which bugs are present in our sample. Now a days, there are many tools and approaches to do that, however there's no gold standard or a consensus database. [Here you can find a review of tools, methods, pipelines and database for analysing microbial communities][3].

In this course, we will use [MetaPhlAn v2.0][4]. In short, MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea, Eukaryotes and Viruses) from metagenomic shotgun sequencing data with strain level resolution. MetaPhlAn 2 relies on ~1M unique specific marker genes identified from ~17,000 reference genomes. Usually, the characterisation of a full metagenomic sample can take ~1h. In order to save some time, we will run a demo. 

1. Check the following [input file][5]. 
2. MetaPhlAn 2 can be run in [Galaxy][6], a web-based platform for data intensive biomedical research. We created a [special session][7] for this course.
3. Open the [link][7] in your browser and click `Import history` in the right upper corner. In the left side of the web page you can see the different available tools, at the right side, the data that we are going to use as a demo, if you click on it you can see the details of the data, download it or edit attributes.  


**Q1. Find `MetaPhlAn 2` in the list of tools in our Galaxy session, read the description and choose the optimal parameters. Execute the program, take a nap and check the output. How many microbial species we are able to identify? What do the numbers represent?**

Import demo data to your session
![Image of first page](https://github.com/ArnauVich/Courses/blob/master/images/page_1_galaxy.png)

![Image Galaxy session](https://github.com/ArnauVich/Courses/blob/master/images/page_2_galaxy.png)

Grey box indicate that the analysis is running. The process can take few minutes, you can refresh your browser to check the status
![Image Galaxy running](https://github.com/ArnauVich/Courses/blob/master/images/page_running.png)

Green!!! It's done! :smile: :smile: 
![Image Galaxy done](https://github.com/ArnauVich/Courses/blob/master/images/page_finished.png)



Module 2: Exploring microbiome data 
---------------------------------------------


Congratulations! You have survived the first module! :muscle: :muscle:

Because we are nice people and you are our favourite students, we decided to provide you a file with 20 microbiome samples already characterized. You can download the files [here][8] and [here][9]

Let's go to R and check these files: 

1. First we will set our working directory and load some packages 
```{r}
	setwd("~/Desktop/Course/")
	library (ggplot2)
	library(scales)
	library(reshape2)
	library (vegan)
	library(foreach)

```


:bangbang::bangbang::bangbang: If you don't have the libraries available you can install them, using: 

```{r}
	install.packages("ggplot")
	
```

2. Import the data

```{r}
	bacteria=read.table("./Microbiome.txt", sep="\t", header=T, row.names = 1)
	phenotypes= read.table("./Phenotypes.txt", sep="\t", header=T, row.names = 1)	
```


**Q2. Check the row names and the column names. What is the structure of each file? Report the number of columns and rows**


3. Next we want to check how many different taxonomies are present in our file (p.e , how many different kingdoms, phyla, species, etc.). A way to do that is using the row names annotation.

```{r}
	> taxas = rownames(bacteria)
	> mini=head(taxas)
	
	#  [TIP!] If we split the now names by "|" and we count the numbers of words in the resulting string we can count the taxomical levels, p.e 1 = Kingdom , 3 = Class

	> strsplit(mini, "\\|")
	
	[[1]]
	[1] "k__Archaea"

	[[2]]
	[1] "k__Archaea"       "p__Euryarchaeota"

	[[3]]
	[1] "k__Archaea"         "p__Euryarchaeota"   "c__Methanobacteria"

	[[4]]
	[1] "k__Archaea"            "p__Euryarchaeota"      "c__Methanobacteria"    "o__Methanobacteriales"

	[[5]]
	[1] "k__Archaea"             "p__Euryarchaeota"       "c__Methanobacteria"     "o__Methanobacteriales"  "f__Methanobacteriaceae"

	[[6]]
	[1] "k__Archaea"             "p__Euryarchaeota"       "c__Methanobacteria"     "o__Methanobacteriales"  "f__Methanobacteriaceae" "g__Methanobrevibacter" 
```


**Q3. How many species, strains, genus, families, etc. we can find in our data? Create a bar plot summarizing the number of different taxonomical levels present in our table** 

***Tip: use a for loop and if conditions***

4. You may also be interested in the mean relative abundance of a microorganism, let's say, how abundant is *Escherichia Coli* in our gut. For that we can simply calculate the mean. In addition we can also calculate only the mean in those samples that the bacteria is present (thus, excluding zeros)

```{r}
	> mean(transposed_bacteria[,1])
	> sum (transposed_bacteria[,1]!=0)
```


**Q4. Create a summary table containing per each bacteria: mean, mean without zero values, and the percentage of samples where it is present. Identify the top 10 most abundant bacteria. How many taxonomies are absent in all the samples**
**Tip: create a matrix to store the results and perform a ***for*** loop**

5. We can use the previous information to remove those microbes that are absent in most of the samples, let's set a threshold of presence of at least 10% of the samples 

```{r}
my_results_filtered=my_results[my_results$perc_missing<90,]
list_to_keep=as.vector(row.names(my_results_filtered))
bacteria_2_keep=bacteria[list_to_keep,]
taxas = rownames(bacteria_2_keep)
```


6. Since the taxonomy is an hierarchical  structure, we may want to perform our analyses only in one specific taxonomical level, let's say species: 


```{r}
list_species=list()
for (i in taxas){
  if (length (unlist(strsplit(i, "\\|"))) == 7){
    list_species=c( list_species,i)
  }
}
species_table=bacteria_2_keep[unlist(list_species), ]
```

**Q5. Create a dataframe containing only phylum level**

7. In the phenotype data frame you can see different information per each sample: age, sex (1:Male, 2:Female), number of sequencing reads, BMI etc. In order to perform a case-control study is better to take into account if there's any difference between groups in other phenotypes that can have an influence on the microbiome composition

**Q6. Plot frequencies or distributions of each phenotype and test if there's any difference between healthy controls and IBD participants**

8. Although tomorrow we are going to perform statistical analyses on the microbiome composition, we can already visualise the differences between groups. Merge taxanomy table created in **Q5** with phenotype table. 

**Q7. Create a stacked bar plot showing the abundance of different phyla comparing different phenotypes (sex, smoking, etc.) and cases vs controls**

9. This plot give us an idea on differences in composition at higher taxonomical levels. However we can also look at interesting bacterial species. 

**Q8. Create boxplots comparing IBD vs Healthy controls of the following species: *Methanobrevibacter smithii*, *Faecalibacterium prausnitzii*, *Escherichia coli* and *Bacteroides vulgatus*. What we can conclude?** 

Example




![Image example rel.abundance species](https://github.com/ArnauVich/Courses/blob/master/images/Rplot.png)

[1]: https://www.lifelines.nl
[2]: https://en.wikipedia.org/wiki/Crohn%27s_disease
[3]: https://github.com/ArnauVich/Courses/blob/master/images/Metagenomics%20-%20Tools%20Methods%20and%20Madness.pdf
[4]: http://huttenhower.sph.harvard.edu/metaphlan2
[5]: https://bitbucket.org/biobakery/humann2/raw/eed75fc7a0d8fe99af8de29ecccea979fc737157/humann2/tests/data/demo.fastq
[6]: https://usegalaxy.org
[7]: http://huttenhower.sph.harvard.edu/galaxy/u/avv/h/coursemetagenomics
[8]: https://github.com/ArnauVich/Courses/blob/master/Microbiome.txt 
[9]: https://github.com/ArnauVich/Courses/blob/master/Phenotypes.txt
