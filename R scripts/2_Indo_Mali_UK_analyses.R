# Introduction
# This part of the study will analyze the taxa differences from the unmapped reads of our Indonesian dataset compared to external global datasets from Mali and the UK. 
# The aim of this analysis is to see what pathogens can be identified in whole blood of sampels in the datasets and if they differ between regions
# We will test this by looking at sample clustering, relative abundance of taxa, differential abundance testing, and diversity estimates. 

## About the samples
# The Indonesian data in this study were generated from the previously-published study by [Natri et al](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008749).
# For this analysis, all RNAseq reads went through quality control using FastQC, then leading and trailing bases below a Phred quality score of 20 and adapters were removed. The remaining reads were then mapped to the human genome (Hg38) using STAR with a two-pass alignment.
# All unmapped reads were then put through further quality control using [KneadData](https://huttenhower.sph.harvard.edu/kneaddata/), where tandem repeats were removed with [TRF](https://tandem.bu.edu/trf/trf.html) and human reads were removed with [BMtagger](https://www.westgrid.ca/support/software/bmtagger).
# Finally, mapping and classification of reads was used with [CCMetagen](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02014-2), which relies upon [KMA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2336-6) for its mapping stage. 

# The control samples in this study are from a study by [Tran et al](https://www.nature.com/articles/srep31291#Tab1), conducted in 2016 and taken from a population in Mali, and by [Singhania et al](https://www.nature.com/articles/s41467-018-04579-w), conducted in 2018 and taken from the UK. 
# The samples in the Tran et al study were whole blood samples from 54 healthy Malian individuals and sequeneced on an Illumina HiSeq 2000. Whole blood was then depleted of rRNA and globin RNA before amplification using the ScriptSeq Complete Gold Kit (Illumina). 
# For the Singhania et al study, whole blood was taken from 10 healthy individuals from London in the United Kingdom. Whole blood was globin-depleted using the human GLOBINclear kit (Thermo Fisher Scientific) and sequenced on an Illumina HiSeq 4000. 
# Both control datasets went through the same pipeline as the Indonesian dataset, however since the Singhania et al data is 75-bp, both the Malian and Indonesian dataset were trimmed to 75BP.
# Therefore, this part of the analysis uses 75-BP, single-ended data from whole blood, taken from three different populations.

# Loading packages, directories, and colour setup ##########
# Install the packages needed to run the analysis.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("pacman")
library("pacman")
install.packages("remotes")
library("remotes")

BiocManager::install(c("phyloseq", "DESeq2", "microbiome", "ALDEx2", "metagenomeSeq", "mixOmics"), lib="/Users/muhamad.fachrul/Library/R/arm64/4.4/library", force=TRUE)
remotes::install_github("malucalle/selbal")
remotes::install_github("microsud/microbiomeutilities")
remotes::install_github("adw96/breakaway")
remotes::install_github("adw96/DivNet")

p_load(RColorBrewer, dplyr, plyr, reshape2, ggpubr, metacoder,
       tidyverse, phyloseq, DESeq2, microbiome, vegan, picante,
       ALDEx2, metagenomeSeq, HMP, dendextend, selbal, rms, breakaway,
       microbiomeutilities, mixOmics, SRS, ggrepel, DivNet, taxize, ggplot2)

# set up directories
refdir = "/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch"

### Make function for aggregate_top_taxa
aggregate_top_taxa <- function (x, top, level) {
  
  x <- aggregate_taxa(x, level)
  
  tops <- top_taxa(x, top)
  tax <- tax_table(x)
  
  inds <- which(!rownames(tax) %in% tops)
  
  tax[inds, level] <- "Other"
  
  tax_table(x) <- tax
  
  tt <- tax_table(x)[, level]
  tax_table(x) <- tax_table(tt)
  
  aggregate_taxa(x, level)
  
}

# set ggplot colour theme to white
theme_set(theme_bw())

# Set up colour schemes
IndonesiaCol="#332288"
MaliCol="#999933"
UKCol="#CC6677"

# Reading in the 75BP Indonesian data ##########
AllREadsPE_Indo_Counts <- read.csv(paste0(refdir,"/Indo_data/CCmetagen/nuRefBM_apm_p_75BP_PE_merged_taxa.csv"),check.names=FALSE)

## Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(AllREadsPE_Indo_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(AllREadsPE_Indo_Counts[,-which(colnames(AllREadsPE_Indo_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

## Convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
AllREadsPE_Indo_Counts_physeq = phyloseq(taxa, tax)

## Add in sample information, starting with Island
samplenames <- colnames(otu_table(AllREadsPE_Indo_Counts_physeq))
pop <- rep("Indonesia",ncol(otu_table(AllREadsPE_Indo_Counts_physeq)))

## Make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(AllREadsPE_Indo_Counts_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(AllREadsPE_Indo_Counts_physeq) <- samples

## Get information on phyloseq objects
AllREadsPE_Indo_Counts_physeq

## We can see that we have a phyloseq object consisting of 123 samples, 6 of which are replicates.
## We'll take the replicates with the highest library depth and then remove the rest.

## Get replicates
trimmed_samplenames = gsub("Batch1",'',samplenames) %>% gsub("Batch2",'', .) %>% gsub("Batch3",'', .)
trimmed_samplenames = sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", trimmed_samplenames)
replicate_index = which(duplicated(trimmed_samplenames) | duplicated(trimmed_samplenames, fromLast = TRUE))
replicates = samplenames[replicate_index]

## Add sequencing depth information to the Physeq object in order to filter replicates by seqDepth
SeqDepth = colSums(otu_table(AllREadsPE_Indo_Counts_physeq))
#sample_data(AllREadsPE_Indo_Counts_physeq)$SeqDepth = SeqDepth

## Find out which replicates have the highest sequencing depth
sample_data(AllREadsPE_Indo_Counts_physeq)[replicates,]
replicateDF=as.data.frame(sample_data(AllREadsPE_Indo_Counts_physeq)[replicates,])
replicateDF$SampleName = sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", replicateDF$SampleName)
replicateDF$SampleName = gsub("Batch1",'',replicateDF$SampleName) %>% gsub("Batch2",'', .) %>% gsub("Batch3",'', .)
replicateDF=replicateDF[with(replicateDF, order(-SeqDepth)), ]
removeReplicates=rownames(replicateDF[which(duplicated(replicateDF$SampleName)),])
keepReplicates=rownames(sample_data(AllREadsPE_Indo_Counts_physeq))[-which(rownames(sample_data(AllREadsPE_Indo_Counts_physeq)) %in% removeReplicates)]

## Prune these out
AllREadsPE_Indo_Counts_physeq=prune_samples(keepReplicates,AllREadsPE_Indo_Counts_physeq)

## Remove taxa with only 0's in the phyloseq object
any(taxa_sums(AllREadsPE_Indo_Counts_physeq) == 0)
AllREadsPE_Indo_Counts_physeq=prune_taxa(taxa_sums(AllREadsPE_Indo_Counts_physeq) > 0, AllREadsPE_Indo_Counts_physeq)

## Add sequencing depth information
SeqDepth_Prefilter = colSums(otu_table(AllREadsPE_Indo_Counts_physeq))
sample_data(AllREadsPE_Indo_Counts_physeq)$SeqDepth_Prefilter = SeqDepth_Prefilter

AllREadsPE_Indo_Counts_physeq

## We now have a phyloseq object of 117 samples

# Reading in the 75BP Malian data ##########
Mali_Counts <- read.csv(paste0(refdir,"/Mali_data/CCmetagen/Mali_nuRefBM_apm_p_75BP_PE_merged_taxa.csv"),check.names=FALSE)

## Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(Mali_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(Mali_Counts[,-which(colnames(Mali_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

## Convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
Mali_Counts_physeq = phyloseq(taxa, tax)

## Add in sample information, i.e., the sample names and population they're from
samplenames <- colnames(otu_table(Mali_Counts_physeq))
pop <- rep("Mali",ncol(otu_table(Mali_Counts_physeq)))

## Make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(Mali_Counts_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(Mali_Counts_physeq) <- samples

## Add sequencing depth information 
SeqDepth_Prefilter = colSums(otu_table(Mali_Counts_physeq))
sample_data(Mali_Counts_physeq)$SeqDepth_Prefilter = SeqDepth_Prefilter

## Summarise data
Mali_Counts_physeq

## We now have a phyloseq object of 54 samples

# Reading in the 75BP UK data ##########
UK_Counts <- read.csv(paste0(refdir,"/UK_data/CCmetagen/UK_nuRefBM_apm_p_75BP_PE_merged_taxa.csv"),check.names=FALSE)

## Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(UK_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(UK_Counts[,-which(colnames(UK_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

## Convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
UK_Counts_physeq = phyloseq(taxa, tax)

## Add in sample information, i.e., the sample names and population they're from
samplenames <- colnames(otu_table(UK_Counts_physeq))
pop <- rep("UK",ncol(otu_table(UK_Counts_physeq)))

## Make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(UK_Counts_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(UK_Counts_physeq) <- samples

## Add sequencing depth information 
SeqDepth_Prefilter = colSums(otu_table(UK_Counts_physeq))
sample_data(UK_Counts_physeq)$SeqDepth_Prefilter = SeqDepth_Prefilter
save(UK_Counts_physeq, file=paste0(refdir,"/UK_data/UK_Counts_physeq.Rda"))

## Get phyloseq summary information
UK_Counts_physeq

## We now have a phyloseq object of 12 samples.

# Merging the data
## Assign unique taxa names to all phyloseq objects
taxa_names(AllREadsPE_Indo_Counts_physeq) <- paste(tax_table(AllREadsPE_Indo_Counts_physeq)[,"Superkingdom"], tax_table(AllREadsPE_Indo_Counts_physeq)[,"Kingdom"], tax_table(AllREadsPE_Indo_Counts_physeq)[,"Phylum"], tax_table(AllREadsPE_Indo_Counts_physeq)[,"Class"], tax_table(AllREadsPE_Indo_Counts_physeq)[,"Order"], tax_table(AllREadsPE_Indo_Counts_physeq)[,"Family"], tax_table(AllREadsPE_Indo_Counts_physeq)[,"Genus"], tax_table(AllREadsPE_Indo_Counts_physeq)[,"Species"], sep="_")
taxa_names(Mali_Counts_physeq) <- make.unique(paste(tax_table(Mali_Counts_physeq)[,"Superkingdom"], tax_table(Mali_Counts_physeq)[,"Kingdom"], tax_table(Mali_Counts_physeq)[,"Phylum"], tax_table(Mali_Counts_physeq)[,"Class"], tax_table(Mali_Counts_physeq)[,"Order"], tax_table(Mali_Counts_physeq)[,"Family"], tax_table(Mali_Counts_physeq)[,"Genus"], tax_table(Mali_Counts_physeq)[,"Species"], sep="_"))
taxa_names(UK_Counts_physeq) <- make.unique(paste(tax_table(UK_Counts_physeq)[,"Superkingdom"], tax_table(UK_Counts_physeq)[,"Kingdom"], tax_table(UK_Counts_physeq)[,"Phylum"], tax_table(UK_Counts_physeq)[,"Class"], tax_table(UK_Counts_physeq)[,"Order"], tax_table(UK_Counts_physeq)[,"Family"], tax_table(UK_Counts_physeq)[,"Genus"], tax_table(UK_Counts_physeq)[,"Species"], sep="_"))
merged_phylo_counts=merge_phyloseq(AllREadsPE_Indo_Counts_physeq, Mali_Counts_physeq, UK_Counts_physeq)

## subset populations
Indonesian_subset <- phyloseq::subset_samples(merged_phylo_counts, SamplePop == "Indonesia")
Mali_subset <- phyloseq::subset_samples(merged_phylo_counts, SamplePop == "Mali")
UK_subset <- phyloseq::subset_samples(merged_phylo_counts, SamplePop == "UK")

## Define function to check if both objects are identical, then check if they are
is.identical <- function(pop_to_subset, original_phyloseq_obj){
  pruned_pop_subset <- prune_taxa(taxa_sums(pop_to_subset) > 0, pop_to_subset)
  merged_df=merge(as.data.frame(otu_table(pruned_pop_subset)), as.data.frame(otu_table(original_phyloseq_obj)),by=0)
  nsamples=ncol(otu_table(pruned_pop_subset))
  merged_df[,"Row.names"] <- NULL
  colnames(merged_df)=rep("Rep",ncol(merged_df))
  is.identical <- identical(merged_df[,1:nsamples],merged_df[,(nsamples+1):ncol(merged_df)])
  return(is.identical)
}

is.identical(pop_to_subset=Indonesian_subset, original_phyloseq_obj=AllREadsPE_Indo_Counts_physeq)
is.identical(pop_to_subset=Mali_subset, original_phyloseq_obj=Mali_Counts_physeq)
is.identical(pop_to_subset=UK_subset, original_phyloseq_obj=UK_Counts_physeq)

# Data processing ##########
## Removing singletons from the data
## Histogram of data
ggplot(meta(merged_phylo_counts)) + geom_histogram(aes(x = SeqDepth_Prefilter), alpha= 0.6, bins=100) + facet_wrap(~SamplePop)

## Separate species' abundances and taxonomy columns
rarecurve_counts <- as.data.frame(otu_table(merged_phylo_counts))
col <- c(rep(IndonesiaCol,sum(sample_data(merged_phylo_counts)[,"SamplePop"]=="Indonesia")),rep(MaliCol,sum(sample_data(merged_phylo_counts)[,"SamplePop"]=="Mali")),rep(UKCol,sum(sample_data(merged_phylo_counts)[,"SamplePop"]=="UK")))
## Try with different filtering thresholds:
for (i in c(1)){
  temp_curvedf <- rarecurve_counts
  temp_curvedf[temp_curvedf<=i]<-0
  rarecurve(t(temp_curvedf), step=200, col=col,label=F, xlab="Counts",ylab="Number Species",main=paste("Removing Reads",i,"And Below",sep=" "),xlim=c(0,50000))
}

temp_curvedf <- rarecurve_counts
temp_curvedf[temp_curvedf<=1]<-0
rarecurve(t(temp_curvedf), step=200, col=col,label=F, xlab="Counts",ylab="Number Species",main=paste("Saturation after singleton removals"),xlim=c(0,50000))

## Rarefaction curves for the Indonesian dataset only (since it is lower sequencing depth and hard to see)
Indonesian_subset <- prune_taxa(taxa_sums(Indonesian_subset) > 0, Indonesian_subset)
temp_Indo_curvedf <- as.data.frame(otu_table(Indonesian_subset))
col <- IndonesiaCol
for (i in c(1)){
  temp_Indo_curvedf[temp_Indo_curvedf<=i]<-0
  rarecurve(t(temp_Indo_curvedf), step=200, col=col,label=F, xlab="Counts",ylab="Number Species",main=paste("Removing Reads",i,"And Below \nIndonesian Samples",sep=" "),xlim=c(0,10000))
}

temp_Indo_curvedf <- as.data.frame(otu_table(Indonesian_subset))
rarecurve(t(temp_Indo_curvedf), step=200, col=IndonesiaCol,label=F, xlab="Counts",ylab="Number of species",main=paste("Indonesia"),xlim=c(0,10000))
temp_UK_curvedf <- as.data.frame(otu_table(UK_subset))
temp_UK_curvedf[temp_UK_curvedf<=1]<-0
rarecurve(t(temp_UK_curvedf), step=200, col=UKCol,label=F, xlab="Counts",ylab="Number of species",main=paste("UK"),xlim=c(0,10000))
rarecurve(t(temp_Indo_curvedf), step=200, col=IndonesiaCol,label=F, xlab="Counts",ylab="Number of species",main=paste("Indonesia"),xlim=c(0,10000)) + rarecurve(t(temp_Mali_curvedf), step=200, col=MaliCol,label=F, xlab="Counts",ylab="Number of species",main=paste("Mali"),xlim=c(0,10000))

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/Indo_curve.png"), units="px", width=2500, height=2000, res=400)
temp_Indo_curvedf[temp_Indo_curvedf<=1]<-0
rarecurve(t(temp_Indo_curvedf), step=200, col=col,label=F, xlab="Counts",ylab="Number Species",main=paste("Indonesia"),xlim=c(0,10000))
#dev.off()

## Filter out singletons
otu_table(merged_phylo_counts)[otu_table(merged_phylo_counts)<=3]<-0
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)
## Add sequencing depth information after filtering out singletons
SeqDepth_noSingletons = colSums(otu_table(merged_phylo_counts))
sample_data(merged_phylo_counts)$SeqDepth_noSingletons = SeqDepth_noSingletons

## Removing humans and plants
### Filter out Viridiplantae 
merged_phylo_counts_noPlant <- subset_taxa(merged_phylo_counts, (Kingdom!="Viridiplantae"))
merged_phylo_counts_noPlant <- prune_taxa(taxa_sums(merged_phylo_counts_noPlant) > 0, merged_phylo_counts_noPlant)
# add sequencing depth information after filtering out plants
SeqDepth_noViridiplantae = colSums(otu_table(merged_phylo_counts_noPlant))
sample_data(merged_phylo_counts_noPlant)$SeqDepth_noViridiplantae = SeqDepth_noViridiplantae

### Filter out Chordata
merged_phylo_counts_noChor <- subset_taxa(merged_phylo_counts_noPlant, (Phylum!="Chordata"))
merged_phylo_counts_noChor <- prune_taxa(taxa_sums(merged_phylo_counts_noChor) > 0, merged_phylo_counts_noChor)
### Add sequencing depth information after filtering out Metazoa
SeqDepth_noChordata = colSums(otu_table(merged_phylo_counts_noChor))
sample_data(merged_phylo_counts_noChor)$SeqDepth_noChordata = SeqDepth_noChordata

### Filter out Metazoa
merged_phylo_counts_noMeta <- subset_taxa(merged_phylo_counts_noChor, (Kingdom!="Metazoa"))
merged_phylo_counts_noMeta <- prune_taxa(taxa_sums(merged_phylo_counts_noMeta) > 0, merged_phylo_counts_noMeta)
### Add sequencing depth information after filtering out Metazoa
SeqDepth_noMetazoa = colSums(otu_table(merged_phylo_counts_noMeta))
sample_data(merged_phylo_counts_noMeta)$SeqDepth_noMetazoa = SeqDepth_noMetazoa

Indonesian_noMz <- phyloseq::subset_samples(merged_phylo_counts_noMeta, SamplePop == "Indonesia")
Indonesian_noMz <- Indonesian_noMz %>%
  tax_glom("Phylum")
Indonesian_noMz <- prune_taxa(taxa_sums(Indonesian_noMz) > 0, Indonesian_noMz)
### 22 taxa, 8 phylum
Indo_ab <- as.data.frame(taxa_sums(Indonesian_noMz)/sum(taxa_sums(Indonesian_noMz)))
Indo_tx <- as.data.frame(taxa_sums(Indonesian_noMz))

Mali_noMz <- phyloseq::subset_samples(merged_phylo_counts_noMeta, SamplePop == "Mali")
Mali_noMz <- Mali_noMz %>%
  tax_glom("Phylum")
Mali_noMz <- prune_taxa(taxa_sums(Mali_noMz) > 0, Mali_noMz)
### 41 taxa, 13 phylum
Mali_ab <- as.data.frame(taxa_sums(Mali_noMz)/sum(taxa_sums(Mali_noMz)))
Mali_tx <- as.data.frame(taxa_sums(Mali_noMz))

UK_noMz <- phyloseq::subset_samples(merged_phylo_counts_noMeta, SamplePop == "UK")
UK_noMz <- UK_noMz %>%
  tax_glom("Phylum")
UK_noMz <- prune_taxa(taxa_sums(UK_noMz) > 0, UK_noMz)
### 101 taxa, 11 phylum
UK_ab <- as.data.frame(taxa_sums(UK_noMz)/sum(taxa_sums(UK_noMz)))
UK_tx <- as.data.frame(taxa_sums(UK_noMz))

### Filter out taxa which are unassigned at the Kingdom level
merged_phylo_counts_noUnk <- subset_taxa(merged_phylo_counts_noMeta, (Superkingdom!="unk_sk"))
merged_phylo_counts_noUnk <- prune_taxa(taxa_sums(merged_phylo_counts_noUnk) > 0, merged_phylo_counts_noUnk)
### Add sequencing depth information after filtering out Metazoa
SeqDepth_noUnkSk = colSums(otu_table(merged_phylo_counts_noUnk))
sample_data(merged_phylo_counts_noUnk)$SeqDepth_noUnkSk = SeqDepth_noUnkSk

All_75BP_PE_postCC <- data.frame(sample_data(merged_phylo_counts_noUnk))
All_75BP_PE_postCC$SampleName <- gsub(".ccm.*","", All_75BP_PE_postCC$SampleName)

#sumdir = "/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/QC/ReadsSummary/"
write.table(All_75BP_PE_postCC[All_75BP_PE_postCC$SamplePop == "Indonesia",], file = paste0(sumdir,"Indo_75BP_PE_postCC.txt"))
write.table(All_75BP_PE_postCC[All_75BP_PE_postCC$SamplePop == "Mali",], file = paste0(sumdir,"Mali_75BP_PE_postCC.txt"))
write.table(All_75BP_PE_postCC[All_75BP_PE_postCC$SamplePop == "UK",], file = paste0(sumdir,"UK_75BP_PE_postCC.txt"))

## Summarising the data
### Barplot of library sizes
#ggplot(meta(merged_phylo_counts_noUnk), aes(SampleName, SeqDepth_noUnkSk)) + geom_bar(stat = "identity", aes(fill = SamplePop)) +
  scale_fill_manual(values = c(IndonesiaCol,MaliCol,UKCol)) + rotate_x_text()
png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/libsizes.png"), units="px", width=2000, height=2000, res=400)
ggplot(meta(merged_phylo_counts_noUnk), aes(SampleName, SeqDepth_noUnkSk)) + geom_bar(stat = "identity", aes(fill = SamplePop)) +
  scale_fill_manual(values = c(IndonesiaCol,MaliCol,UKCol)) + rotate_x_text()
#dev.off()

## We can see that the library sizes are highly uneven, with the Indonesian data having the lowest sampling depth (with the exception of a few samples) and the UK dataset having the highest library depth.
## The final step us is to summarise the data and see how many reads we lost at each filtering step.
FilteringSummary <- sample_data(merged_phylo_counts_noUnk)[,c("SamplePop","SeqDepth_Prefilter","SeqDepth_noSingletons","SeqDepth_noViridiplantae","SeqDepth_noChordata","SeqDepth_noMetazoa", "SeqDepth_noUnkSk")]

### Melt df and plot
melted_FilteringSummary = melt(FilteringSummary)
ggplot(melted_FilteringSummary, aes(x=variable, y=value, fill=SamplePop)) +
  geom_violin(alpha=0.8) + theme_bw() + ylab("Sampling depth") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = c(IndonesiaCol,MaliCol,UKCol)) +
  geom_boxplot(color="black",width=0.2, alpha = 0.7) + facet_wrap(~ SamplePop, scales = "free")

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/sampDepth.png"), units="px", width=2500, height=2000, res=400)
ggplot(melted_FilteringSummary, aes(x=variable, y=value, fill=SamplePop)) +
  geom_violin(alpha=0.8) + theme_bw() + ylab("Sampling depth") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = c(IndonesiaCol,MaliCol,UKCol)) +
  geom_boxplot(color="black",width=0.2, alpha = 0.7) + facet_wrap(~ SamplePop, scales = "free")
#dev.off()

# Data normalisation ##########
## Centered log-ration transformation
## Taxa can be viewed by their relative abundance, however changes in the abundance of one taxon will result in changing the abundance of other taxa.
## One of the ways to handle this is to transform the data using Centered Log Ratio (CLR) transformation. CLR data shows how OTUs behave relative to the per-sample average and is a commonly-used data transformation method in microbiomics.

## CLR transformation is quite sensitive to zeros, so we need to remove as many zeros as we can before performing CLR transformation.
## We'll do this by merging our taxa data at the phylum level, and add an offset of 0.1.
pop_comparison <- merged_phylo_counts_noUnk %>%
  tax_glom("Phylum")

pop_comparison

## Change taxa names to reflect Phylum level
taxa_names(pop_comparison) <- paste(tax_table(pop_comparison)[,"Superkingdom"], tax_table(pop_comparison)[,"Kingdom"], tax_table(pop_comparison)[,"Phylum"], sep="_")

zComposition_estimate <- otu_table(pop_comparison) + 1
zComposition_clr <- microbiome::transform(zComposition_estimate, "clr")

## Add in zCompositions information to new phyloseq object
merged_phylo_counts_zComposition <- pop_comparison
taxa_zComposition <- otu_table(zComposition_clr, taxa_are_rows = TRUE)
otu_table(merged_phylo_counts_zComposition) <- taxa_zComposition

## Sample grouping
### Make an ordination plot using euclidean distances
zComposition.ord <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")

### plot by Island
pcs_vase <- cbind(zComposition.ord$vectors, sample_data(merged_phylo_counts_zComposition))

cohorts_PCA1_2 <- ggplot(pcs_vase) + geom_point(aes(x=Axis.1, y=Axis.2, color=SamplePop, shape=SamplePop, ), alpha=0.6, size=2.5) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_shape_manual(name = "Population", values = c(16, 15, 17)) + scale_colour_manual(name = "Population", values=c(IndonesiaCol,MaliCol,UKCol)) + labs(color = "Population")
cohorts_PCA1_2

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/cohorts_PCs1_2.png"), units="px", width=2500, height=2000, res=400)
cohorts_PCA1_2
#dev.off()

plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 2:3, label="SampleName") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))

cohorts_PCA3_4 <- plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 3:4) + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + labs(color = "Population")
cohorts_PCA3_4

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/cohorts_PCs3_4.png"), units="px", width=2500, height=2000, res=400)
cohorts_PCA3_4
#dev.off()

plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 4:5, label="SampleName") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="SamplePop", axes = 5:6, label="SampleName") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))

## Plot by Plasmodiidae
logged_phyla_counts = log10(colSums(otu_table(merged_phylo_counts_noUnk)[grep("Plasmodiidae",tax_table(merged_phylo_counts_noUnk)[,"Family"])])+1)
sample_data(merged_phylo_counts_zComposition)[["Plasmodiidae"]] = logged_phyla_counts
out.wuf.log <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")

pcs_vase_plas <- cbind(zComposition.ord$vectors, sample_data(merged_phylo_counts_zComposition))

Plas_PCs1_2 <- plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Plasmodiidae", axes = 1:2, label="SampleName") + geom_point(aes(shape = SamplePop), alpha=0.6, size=4) + scale_shape_manual(values = c(16, 15, 17)) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#A03F03"))(10), limits=c(1, 6))
#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/Api_PCs1_2.png"), units="px", width=2500, height=2000, res=400)
Plas_PCs1_2
#dev.off()

Plas_PCs3_4 <- ggplot(pcs_vase_plas) + geom_point(aes(x=Axis.3, y=Axis.4, color=Plasmodiidae, shape=SamplePop), alpha=0.7, size=2.5) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_shape_manual(name = "Population", values = c(16, 15, 17), guide = "none") + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#A03F03"))(10), limits=c(1, 6))
#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/Api_PCs3_4.png"), units="px", width=2500, height=2000, res=400)
Plas_PCs3_4
#dev.off()

ggplot(pcs_vase_plas) + geom_point(aes(x=Axis.5, y=Axis.6, color=Plasmodiidae, shape=SamplePop), alpha=0.7, size=2.5) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_shape_manual(name = "Population", values = c(16, 15, 17), guide = "none") + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#A03F03"))(10), limits=c(1, 6))

cor(pcs_vase_plas[,1:2], log10_plasm_counts, method = "pearson") #Correlation between PCs and load

## Plot by Flaviviridae
logged_phyla_counts = log10(colSums(otu_table(merged_phylo_counts_noUnk)[grep("Flaviviridae",tax_table(merged_phylo_counts_noUnk)[,"Family"])])+1)
sample_data(merged_phylo_counts_zComposition)[["Flaviviridae"]] = logged_phyla_counts
out.wuf.log <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")

pcs_vase_flav <- cbind(zComposition.ord$vectors, sample_data(merged_phylo_counts_zComposition))

Flav_PCs1_2 <- plot_ordination(merged_phylo_counts_zComposition, zComposition.ord, color="Flaviviridae", axes = 1:2, label="SampleName") + geom_point(aes(shape = SamplePop), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/Flav_PCs1_2.png"), units="px", width=2500, height=2000, res=400)
Flav_PCs1_2
#dev.off()

Flav_PCs3_4 <- ggplot(pcs_vase_flav) + geom_point(aes(x=Axis.3, y=Axis.4, color=Flaviviridae, shape=SamplePop, ), alpha=0.7, size=2.5) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_shape_manual(name = "Population", values = c(16, 15, 17), guide = "none") + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))
#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/Flav_PCs3_4.png"), units="px", width=2500, height=2000, res=400)
Flav_PCs3_4
#dev.off()

ggplot(pcs_vase_flav) + geom_point(aes(x=Axis.5, y=Axis.6, color=Flaviviridae, shape=SamplePop, ), alpha=0.7, size=2.5) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_shape_manual(name = "Population", values = c(16, 15, 17), guide = "none") + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 4))

cor(pcs_vase_flav[,1:2], log10_flav_counts, method = "pearson") #Correlation between PCs and load

## Hierarchical clustering by Euclidean distance
ps_otu <- data.frame(phyloseq::otu_table(merged_phylo_counts_zComposition))
ps_otu <- t(ps_otu)
bc_dist <- vegan::vegdist(ps_otu, method = "euclidean")
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
### Provide color codes
meta <- data.frame(phyloseq::sample_data(merged_phylo_counts_zComposition))
colorCode <- c(Indonesia = IndonesiaCol, `Mali` = MaliCol, `UK` = UKCol)
labels_colors(ward) <- colorCode[meta$SamplePop][order.dendrogram(ward)]
### Plot
#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/cohorts_phylo.png"), units="px", width=1500, height=1500, res=320)
par(mar = c(5, 5, 5, 5), cex=0.4)
plot(ward)
#dev.off()

# Relative frequency of taxa  ##########
## Add a new column containing family names and superkingdom
tax_table(merged_phylo_counts_noUnk)[,"Superkingdom"] <- paste(tax_table(merged_phylo_counts_noUnk)[,"Superkingdom"], tax_table(merged_phylo_counts_noUnk)[,"Family"], sep="_")

tax_table(merged_phylo_counts_noUnk)[,"Superkingdom"] <- str_replace_all(tax_table(merged_phylo_counts_noUnk)[,"Superkingdom"], "Bacteria_unk_f", "Bacteria_unclassified")

tax_table(merged_phylo_counts_noUnk)[,"Superkingdom"] <- str_replace_all(tax_table(merged_phylo_counts_noUnk)[,"Superkingdom"], "Eukaryota_unk_f", "Eukaryota_unclassified")

tax_table(merged_phylo_counts_noUnk)[,"Superkingdom"] <- str_replace_all(tax_table(merged_phylo_counts_noUnk)[,"Superkingdom"], "Viruses_unk_f", "Viruses_unclassified")

aggregated_phyloCounts <- aggregate_top_taxa(merged_phylo_counts_noUnk, "Superkingdom", top = 20)
# transform to relative counts
relative_phyloCounts <- microbiome::transform(aggregated_phyloCounts, "compositional")
# Remove weird extra family names added at the end of Superkingdom names
tax_table(relative_phyloCounts)[,"Superkingdom"] <- paste(sapply(strsplit(taxa_names(relative_phyloCounts), "[_.]"), `[`, 1), sapply(strsplit(taxa_names(relative_phyloCounts), "[_.]"), `[`, 2), sep="_")
# Change "Other_NA" to just "Other"
tax_table(relative_phyloCounts)[,"Superkingdom"][grep("Other", taxa_names(relative_phyloCounts))] = "Other"


## Plot
p=plot_bar(relative_phyloCounts, fill = "Superkingdom")

## Set colour palette
families=levels(as.factor(p$data$Superkingdom))
## Get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(10)
PaletteEukaryote = colorRampPalette(c("#4D1C00", "#FF8E4D"))(6)
#PaletteEukaryote <- c(colorRampPalette(c("#b20000","#ff4c4c"))(6)[1:2],"gold",colorRampPalette(c("#b20000","#ff4c4c"))(6)[4:6])
#"#800026","#fd8d3c"
PaletteOther = colorRampPalette(c("black"))(1)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(4)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteOther,PaletteVirus)

phyloseq::plot_bar(relative_phyloCounts, fill = "Superkingdom") +
  geom_bar(aes(fill = Superkingdom), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") + theme_bw(base_size = 15) +
  facet_wrap(~ SamplePop, scales = "free") + scale_fill_manual(values=Merged_Palette, na.translate = F) + theme(axis.text.x = element_text(angle = 90))
theme(panel.background = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/cohorts_relative_ab_gold.png"), units="px", width=3500, height=2000, res=320)
phyloseq::plot_bar(relative_phyloCounts, fill = "Superkingdom") +
  geom_bar(aes(fill = Superkingdom), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") + theme_bw(base_size = 15) +
  facet_wrap(~ SamplePop, scales = "free") + scale_fill_manual(values=Merged_Palette, na.translate = F) + theme(axis.text.x = element_text(angle = 90)) + theme(panel.background = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
#dev.off()

## Let's plot this using a circular barplot
data=as.matrix(as.data.frame(otu_table(relative_phyloCounts)))
data=t(data)
data=as.data.frame(data)
data$group=as.character(as.matrix(as.data.frame(phyloseq::sample_data(merged_phylo_counts_noUnk)[,"SamplePop"])))
data$individual=rownames(data)
data$individual=as.factor(data$individual)
data$group=as.factor(data$group)

tax_rep <- tax_table(relative_phyloCounts)[,"Superkingdom"]

colnames(data)[1:21] <- tax_rep

## Transform data in a tidy format (long format)
data = data %>% gather(key = "observation", value="value", -c(22,23)) 

## Set a number of 'empty bar' to add at the end of each group
empty_bar=6
nObsType=nlevels(as.factor(data$observation))
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar*nObsType )
data=rbind(data, to_add)
data=data %>% arrange(group, individual)
data$id=rep( seq(1, nrow(data)/nObsType) , each=nObsType)

## Prepare a data frame for base lines
base_data=data %>% 
  group_by(group) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

## Prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1

## Make the plot
p = ggplot(data) +      
  
  ### Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.8) +
  # Add a valu=100/75/50/25 lines. I do it at the beginning to make sure barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.25, xend = start, yend = 0.25), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.50, xend = start, yend = 0.50), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.75, xend = start, yend = 0.75), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1, xend = start, yend = 1), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  ### Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),5), y = c(0, 0.25, 0.50, 0.75, 1), label = c("0", "0.25", "0.50", "0.75", "1") , color="black", size=2 , angle=0, fontface="bold", hjust=1) +
  
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=base_data, aes(x = title, y = 1.2, label=group), hjust=c(0,0,0), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE) +
  scale_fill_manual(values=Merged_Palette, na.translate = F)

p

patho_plt <- (Plas_PCs3_4 + Flav_PCs3_4) + plot_layout(guides = 'collect', nrow = 1) & theme(legend.box = "horizontal", legend.text = element_text(size=6), legend.key.size = unit(0.5, 'cm'))

p / patho_plt + plot_layout(heights = c(1.8, 1))

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/circularGlobal.png"), units="px", width=2500, height=1500, res=300)
pdf(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/circularGlobal.pdf"), width=8, height=4)
p + theme(legend.box = "vertical", legend.text = element_text(size=6), legend.key.size = unit(0.4, 'cm'), legend.title = element_text(size=6))
dev.off()

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/cohorts_fig1_gold.png"), units="px", width=2500, height=2000, res=400)
p / patho_plt + plot_layout(heights = c(1.8, 1)) & theme(legend.box = "horizontal", legend.text = element_text(size=6), legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=8))
#dev.off()

# Differential abundance testing ##########
## We're interested in testing whether the taxa composition between populations is significantly different.
## Aldex2 corrects for uneven library sizes, so rarefying is not necessary. The only input needed is the data with singletons removed.
## Since Aldex2 works pairwise, the comparison is done first between the Indonesian and the UK datasets, then the Indonesian and the Malian datasets.

## UK vs Indo
UKvIndo <- subset_samples(pop_comparison, SamplePop != "Mali")
### Remove any 0s from the data
UKvIndo <- prune_taxa(taxa_sums(UKvIndo) > 0, UKvIndo)
taxa_names(UKvIndo)=make.unique(tax_table(UKvIndo)[,"Phylum"])

UKvIndo <- prune_samples(sample_sums(UKvIndo) > 0, UKvIndo)

### Run aldex2
aldex2_UKvIndo <- ALDEx2::aldex(data.frame(phyloseq::otu_table(UKvIndo)), phyloseq::sample_data(UKvIndo)$SamplePop, test="t", effect = TRUE)

### Let's now pull out the significant values after getting the ALdex2 ouput. We're interested in the columns 
sig_aldex2_UKvIndo <- aldex2_UKvIndo %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

### Add in taxonomic information to the significant UKvIndo object
taxa_info <- data.frame(tax_table(UKvIndo)[,c("Superkingdom","Kingdom","Phylum")])
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
sig_aldex2_UKvIndo <- left_join(sig_aldex2_UKvIndo, taxa_info)

### Now for the non-significant df
aldex2_UKvIndo$OTU <- rownames(aldex2_UKvIndo)
aldex2_UKvIndo <- left_join(aldex2_UKvIndo, taxa_info)

### Set significance colours
#aldex2_UKvIndo$threshold <- aldex2_UKvIndo$we.eBH <= 0.05
#aldex2_UKvIndo$threshold = as.numeric(aldex2_UKvIndo$threshold) + 1

### Add label names
labels_uvi = c()
counter = 0
for (name in as.character(aldex2_UKvIndo$Phylum)){
  counter = counter + 1
  if (name == "unk_p"){
    name = as.character(aldex2_UKvIndo$Superkingdom[counter])
  }
  labels_uvi = c(labels_uvi, name)
}
### Get superkingdom names
taxa_superkingdom = sapply(strsplit(tax_table(UKvIndo)[,"Superkingdom"], "[_.]"), `[`, 1)

### Plot
UKvIndo_plot <- ggplot(aldex2_UKvIndo) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_UKvIndo$we.eBH <= 0.05, c("#023858","#800026","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=2.5) +
  #geom_text_repel(aes(x = effect, y = -log10(wi.eBH), label = rownames(aldex2_IndoVsDutch))) +
  ggtitle("UK versus Indonesia") +
  xlab("Effect Size") + 
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  geom_vline(xintercept=c(-0.6, 0.6), col="darkgrey") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.05, labels_uvi,""))) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/UK_vs_Indo_DA.png"), units="px", width=2500, height=2000, res=400)
UKvIndo_plot
#dev.off()

## UK vs Indo
MalivIndo=subset_samples(pop_comparison, SamplePop != "UK")
any(taxa_sums(MalivIndo) == 0)
MalivIndo <- prune_taxa(taxa_sums(MalivIndo) > 0, MalivIndo)
taxa_names(MalivIndo)=make.unique(tax_table(MalivIndo)[,"Phylum"])

MalivIndo <- prune_samples(sample_sums(MalivIndo) > 0, MalivIndo)

aldex2_MalivIndo <- ALDEx2::aldex(data.frame(phyloseq::otu_table(MalivIndo)), phyloseq::sample_data(MalivIndo)$SamplePop, test="t", effect = TRUE)
sig_aldex2_MalivIndo <- aldex2_MalivIndo %>%
  rownames_to_column(var = "OTU") %>%
  filter(we.eBH < 0.05) %>%
  arrange(effect, we.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)

# add in taxonomic information to the significant MalivIndo object
taxa_info <- data.frame(tax_table(MalivIndo)[,c("Superkingdom","Kingdom","Phylum")])
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
sig_aldex2_MalivIndo <- left_join(sig_aldex2_MalivIndo, taxa_info)
# now for the non-significant df
aldex2_MalivIndo$OTU <- rownames(aldex2_MalivIndo)
aldex2_MalivIndo <- left_join(aldex2_MalivIndo, taxa_info)

# set significance colours
aldex2_MalivIndo$threshold <- aldex2_MalivIndo$we.eBH <= 0.05
aldex2_MalivIndo$threshold = as.numeric(aldex2_MalivIndo$threshold) + 1

# add label names
labels_mvi = c()
counter = 0
for (name in as.character(aldex2_MalivIndo$Phylum)){
  counter = counter + 1
  if (name == "unk_p"){
    name = as.character(aldex2_MalivIndo$Superkingdom[counter])
  }
  labels_mvi = c(labels_mvi, name)
}
# get superkingdom names
taxa_superkingdom = sapply(strsplit(tax_table(MalivIndo)[,"Superkingdom"], "[_.]"), `[`, 1)

# plot
MalivIndo_plot <- ggplot(aldex2_MalivIndo) +
  geom_point(aes(x = effect, y = -log10(we.eBH)), color = ifelse(aldex2_MalivIndo$we.eBH <= 0.05, c("#023858","#800026","#78c679")[as.numeric(as.factor(taxa_superkingdom))],"black"), alpha = 0.65, size=2.5) +
  ggtitle("Mali versus Indonesia") +
  xlab("Effect Size") +
  ylab("-log10 adjusted p-value") + theme_bw(base_size = 18) +
  geom_vline(xintercept=c(-0.6, 0.6), col="darkgrey") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  geom_text_repel(aes(x = effect, y = -log10(we.eBH), label = ifelse(we.eBH <= 0.05, labels_mvi,""))) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/1.PE_IndonesianVsMaliAndTBControls_CLRMethod_KeepingOutlier_files/Mali_vs_Indo_DA.png"), units="px", width=2500, height=2000, res=400)
MalivIndo_plot
#dev.off()
