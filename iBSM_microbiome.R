setwd("~/Desktop/Course/")
library (ggplot2)
library(scales)
library(reshape2)
library (vegan)
library(foreach)
#Upload data
bacteria=read.table("./Microbiome.txt", sep="\t", header=T, row.names = 1)
phenotypes= read.table("./Phenotypes.txt", sep="\t", header=T, row.names = 1)
# Trick to check the data: Samples in columns, taxonomy in rows 
head (row.names(bacteria))
head (colnames(bacteria))
dim (bacteria)
# Set variables to zero
kingdoms=0
phylum=0
class=0
order=0
family=0
genus=0
species=0
strains=0
# List all rownames
taxas = rownames(bacteria)
# Get first 5
mini=head(taxas)
# Split pipes 
strsplit(mini, "\\|")
# Count different levels. 
for (i in taxas){
  if (length (unlist(strsplit(i, "\\|"))) == 1){
    kingdoms=kingdoms +1
  }else if (length (unlist(strsplit(i, "\\|"))) == 2) {
    phylum=phylum +1
  } else if (length (unlist(strsplit(i, "\\|"))) == 3) {
    class=class +1
  } else if (length (unlist(strsplit(i, "\\|"))) == 4) {
    order=order +1
  } else if (length (unlist(strsplit(i, "\\|"))) == 5) {
    family=family +1
  } else if (length (unlist(strsplit(i, "\\|"))) == 6) {
    genus=genus +1 
  } else if (length (unlist(strsplit(i, "\\|"))) == 7) {
    species=species +1
  } else if (length (unlist(strsplit(i, "\\|"))) == 8){
    strains=strains +1
  }
}
#Put results in a table
summary_taxa=as.data.frame(c(kingdoms,phylum,class,order,family,genus,species,strains), row.names = c("1_kingdoms","2_phylum","3_class","4_order","5_family","6_genus","7_species","8_strains"))
colnames(summary_taxa)="counts"
#Plot counts
ggplot(summary_taxa, aes(rownames(summary_taxa),counts)) + geom_bar(stat = "identity") +theme_classic()

##Traspose to calculate stats per columns 
transposed_bacteria=as.data.frame(t(bacteria))

##Let's calculate the mean values per the first taxa, in this case, Kingdom Archaea
mean(transposed_bacteria[,1])

## Check the number of non 0's values
sum (transposed_bacteria[,1]!=0)

my_results=matrix(ncol = 5, nrow=ncol(transposed_bacteria)) 
taxonomy_abundance <- function(taxonomy_table) {
  ##Function to calculate mean excluding 0 values
  nzmean <- function(a){
    mean(a[a!=0])
  }
  ##Function to calculate nº of 0
  zsum <- function(a){
    sum (a==0)
  }
  ##Function to calculate nº of non-0
  nsum <- function(a){
    sum (a!=0)
  }
  ## Loop for each column (taxonomy) in the taxonomy table
  for (i in 1:ncol(taxonomy_table)) {
    #Calculate mean for each column
    aa = mean(taxonomy_table[,i])
    #Calculate number of non-zeros (individuals)
    bb = nsum(taxonomy_table[,i])
    #Calculate mean without taking into account the 0
    cc = nzmean(taxonomy_table[,i])
    #Calculate number of zeros 
    dd = zsum(taxonomy_table[,i])
    ee= (dd/(dd+bb))*100
    my_results[i,1] = aa
    my_results[i,2] = bb
    my_results[i,3] = cc
    my_results[i,4] = dd
    my_results[i,5] = ee
  }
  return(my_results)
}
my_results=as.data.frame(taxonomy_abundance(transposed_bacteria))
rownames(my_results) = colnames(transposed_bacteria)
colnames(my_results) = c("Mean","N_of_non-0", "Non-0_Mean", "N_of_0", "perc_missing") 

# Plot percentage of missingness
ggplot(my_results, aes(perc_missing)) + geom_bar() +theme_classic()

my_results_filtered=my_results[my_results$perc_missing<90,]
list_to_keep=as.vector(row.names(my_results_filtered))
bacteria_2_keep=bacteria[list_to_keep,]
taxas = rownames(bacteria_2_keep)
#Filter 
list_species=list()
for (i in taxas){
  if (length (unlist(strsplit(i, "\\|"))) == 7){
    list_species=c( list_species,i)
  }
}
species_table=bacteria_2_keep[unlist(list_species), ]

list_phylum=list()
for (i in taxas){
  if (length (unlist(strsplit(i, "\\|"))) == 2){
    list_phylum=c( list_phylum,i)
  }
}
phyla_table=bacteria_2_keep[unlist(list_phylum), ]
#phyla_results=my_results[unlist(list_phylum), ]

phyla_pheno=merge(phenotypes, t(phyla_table), by="row.names")
phyla_pheno$Sex=NULL
phyla_pheno$PFReads=NULL
phyla_pheno$Age=NULL
phyla_pheno$BMI=NULL
phyla_pheno$Smoking=NULL
phyla_pheno$PPI=NULL
melt_phyla=melt(phyla_pheno)
ggplot(melt_phyla, aes(DiagnosisCurrent, as.numeric (value))) + geom_bar (aes(fill = variable), stat = "identity") + theme_classic() + xlab("Group") + ylab("relative_abundance")  + scale_y_continuous(labels = percent_format())

#Calculate Alpha-diversity Shannon index (go back to sp)

transposed_species=as.data.frame(t(species_table))
alpha_shannon <- as.data.frame(diversity(transposed_species, index="shannon"))
colnames(alpha_shannon)[1] <- "Shannon_Index"
num_sp=data.frame(SID=character(), Number_of_Species=integer(), stringsAsFactors=FALSE)
for (i in 1:nrow(transposed_species)) {
  num_sp[i,1]=rownames(transposed_species)[i]
  num_sp[i,2]=sum(transposed_species[i,]!=0)
}
num_sp=as.data.frame(num_sp)
row.names(num_sp)=num_sp$SID
colnames(num_sp)=c("SID","Number_of_species")
phenotypes_with_diversity=merge(phenotypes,num_sp, by="row.names")
rownames(phenotypes_with_diversity)=phenotypes_with_diversity$Row.names
phenotypes_with_diversity$Row.names=NULL

phenotypes_with_diversity=merge(phenotypes_with_diversity,alpha_shannon, by="row.names")
rownames(phenotypes_with_diversity)=phenotypes_with_diversity$Row.names
phenotypes_with_diversity$Row.names=NULL
phenotypes_with_diversity$SID=NULL

ggplot(phenotypes_with_diversity, aes(DiagnosisCurrent, as.numeric(Number_of_species), fill = DiagnosisCurrent)) + labs (y="Number_of_species", x="Category") + geom_violin(trim=FALSE) + geom_boxplot(width = 0.1) + theme_classic() + theme(legend.position="none") + theme(axis.text.x = element_text(hjust = 1, size=16,color="black")) + ylim(0,150)

ggplot(phenotypes_with_diversity, aes(DiagnosisCurrent, as.numeric(Shannon_Index), fill = DiagnosisCurrent)) + labs (y="Shannon_Index", x="Category") + geom_violin(trim=FALSE) + geom_boxplot(width = 0.1) + theme_classic() + theme(legend.position="none") + theme(axis.text.x = element_text(hjust = 1, size=16,color="black"))

wilcox.test(phenotypes_with_diversity$Number_of_species ~ phenotypes_with_diversity$DiagnosisCurrent)
wilcox.test(phenotypes_with_diversity$Shannon_Index ~ phenotypes_with_diversity$DiagnosisCurrent)

## Check correlation between phenotypes and bacterial divesity
cor.test (phenotypes_with_diversity$Age, phenotypes_with_diversity$Number_of_species)
plot (phenotypes_with_diversity$Age, phenotypes_with_diversity$Number_of_species)
## Are the results the same when we test the correlations per each cohort?
IBD=subset(phenotypes_with_diversity,phenotypes_with_diversity$DiagnosisCurrent=="IBD")
Controls=subset(phenotypes_with_diversity,phenotypes_with_diversity$DiagnosisCurrent=="Controls")
cor.test(IBD$Age, IBD$Number_of_species, method = "spearman")
#beta diversity

beta<- vegdist(transposed_species, method="bray")
mypcoa=cmdscale(beta, k = 5)
mypcoa2=as.data.frame(mypcoa)
plot(mypcoa2$V1,mypcoa2$V2)

phenotypes_with_diversity=merge(phenotypes_with_diversity,mypcoa2,by="row.names")
phenotypes_with_diversity$color="black"
phenotypes_with_diversity[  phenotypes_with_diversity$DiagnosisCurrent=="IBD",]$color <- "red"
ggplot (phenotypes_with_diversity, aes(V1, V2, geom="blank", colour=color)) + geom_point () + scale_color_identity ("Datasets", breaks=phenotypes_with_diversity$color, labels=phenotypes_with_diversity$DiagnosisCurrent, guide="legend") + theme_classic() + labs(x="PCoA1", y="PCoA2")

example_NMDS=metaMDS(transposed_species, k=2)
stressplot(example_NMDS)
plot(example_NMDS)
ordiplot(example_NMDS,type="n")
##orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)
ordihull(example_NMDS,groups=phenotypes_with_diversity$DiagnosisCurrent,draw="polygon",col="grey90",label=F)
fit <- envfit(example_NMDS, phenotypes , perm = 999) 
plot(example_NMDS)
plot(fit)

example_NMDS=metaMDS(t(phyla_table), k=2)
stressplot(example_NMDS)
plot(example_NMDS)
ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)
ordihull(example_NMDS,groups=phenotypes_with_diversity$DiagnosisCurrent,draw="polygon",col="grey90",label=F)

## Adonis (Diagnosis and Number of seq. reads p-value significant)
adon<-foreach(i=1:ncol(phenotypes),.combine=rbind)%do%{
    
         ad1<-adonis(beta ~ phenotypes[,i],permutations=999)
         ad1$aov.tab[1,]
     }
rownames(adon) = colnames(phenotypes)

## Maybe add cladogram? Ideas? 

## Explore the distribution of the rel. abundance random bacteria

# We see 1) Large number of zeros 2) No normally distributed data

hist(transposed_bacteria[,200], breaks=56)
hist(transposed_bacteria[,2], breaks=60)

#Info on data transformation 
#http://www.biostathandbook.com/transformation.html
#http://strata.uga.edu/8370/rtips/proportions.html

# First we devided values by 100 to adjust scale 0-1 

divided_sp=transposed_species/100

