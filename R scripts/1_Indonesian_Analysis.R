# Introduction
# This part of the study will analyse differences unmapped reads in our Indonesian samples. 
# The aim of this analysis is therefore to see what pathogen can be identified in whole blood in our Indonesian samples and if it differs between islands.
# We will test this by looking at sample clustering, relative abundance of taxa, differential abundance testing, and diversity estimates. 

## About the samples
# The Indonesian data in this study were generated from the previously-published study by [Natri et al](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008749).
# For this analysis, all RNAseq reads went through quality control using FastQC, then leading and trailing bases below a Phred quality score of 20 and adapters were removed. The remaining reads were then mapped to the human genome (Hg38) using STAR with a two-pass alignment.
# All unmapped reads were then put through further quality control using [KneadData](https://huttenhower.sph.harvard.edu/kneaddata/), where tandem repeats were removed with [TRF](https://tandem.bu.edu/trf/trf.html) and human reads were removed with [BMtagger](https://www.westgrid.ca/support/software/bmtagger).
# Finally, mapping and classification of reads was used with [CCMetagen](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02014-2), which relies upon [KMA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2336-6) for its mapping stage. 


# Loading packages, directories, and colour setup ##########
# Install the packages needed to run the analysis.
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(reshape2)
library(ggpubr)
library(metacoder)
library(tidyverse)             
library(phyloseq)                   
library(DESeq2)                       
library(microbiome)               
library(vegan)                         
library(picante)                     
library(ALDEx2)                      
library(metagenomeSeq)          
library(HMP)                             
library(dendextend)               
library(selbal)                       
library(rms)
library(breakaway)        
library(microbiomeutilities)
library(mixOmics)
library(SRS)
library(ggrepel)
library(DivNet)
library(zCompositions)
library(factoextra)
library(compositions)
library(patchwork)
library(edgeR)

## Set up directories
refdir = "/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/CCmetagen/"
batchInfodir = "/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/EpiStudy/forBatch/"

## Make function for aggregate_top_taxa
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

## Set ggplot colour theme to white
theme_set(theme_bw())

## Set up colour schemes
KorowaiCol="#F21A00"
MentawaiCol="#3B9AB2"
SumbaCol="#EBCC2A"

# Reading in the 101BP Indonesian data ##########
AllREadsPE_Indo_Counts <- read.csv(paste0(refdir,"nuRefBM_apm_p_101BP_PE_merged_taxa.csv"),check.names=FALSE)

## Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(AllREadsPE_Indo_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(AllREadsPE_Indo_Counts[,-which(colnames(AllREadsPE_Indo_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

## Convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
AllREadsPE_Indo_Counts_physeq = phyloseq(taxa, tax)

## Add in sample information, starting with Island
samplenames <- colnames(otu_table(AllREadsPE_Indo_Counts_physeq))
pop <- sapply(strsplit(samplenames, "[-.]"), `[`, 1)
pop <- str_replace_all(pop, "MPI", "KOR")

## Add in batch information
load(paste0(batchInfodir, "dataPreprocessingindoRNA.read_counts.TMM.filtered.Rda"))
batch_df = data.frame(rownames(y$samples),y$samples$batch)
colnames(batch_df)=c("Sample","Batch")
batch_df=batch_df[order(batch_df$Sample),]
batch=batch_df[,"Batch"]

## Make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(AllREadsPE_Indo_Counts_physeq)), SamplePop=pop, batch=batch)
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
sample_data(AllREadsPE_Indo_Counts_physeq)$SeqDepth = SeqDepth

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
cat("After removing replicates:\n")
AllREadsPE_Indo_Counts_physeq

## We now have a phyloseq object of 117 samples and 1,172 taxa.

# Data processing ##########
## Removing singletons from the data
# histogram of data
ggplot(meta(AllREadsPE_Indo_Counts_physeq)) + geom_histogram(aes(x = SeqDepth_Prefilter), alpha= 0.6, bins=100)

## The sequencing depth is variable and quite low in some samples.
## Therefore, pushing this threshold up too high will eliminate rare taxa, especially given that we didn't have a high library size to begin with.
## We will stick with removing singletons. We will also add the sequencing depth information to the phyloseq object to keep track of the library size after filtering.

## Save unfiltered phyloseq object
AllREadsPE_Indo_Counts_physeq_withSingletons <- AllREadsPE_Indo_Counts_physeq
## Filter out singletons
AllREadsPE_Indo_Counts_physeq_noSing <- AllREadsPE_Indo_Counts_physeq
otu_table(AllREadsPE_Indo_Counts_physeq_noSing)[otu_table(AllREadsPE_Indo_Counts_physeq_noSing)<=1]<-0
AllREadsPE_Indo_Counts_physeq_noSing <- prune_taxa(taxa_sums(AllREadsPE_Indo_Counts_physeq_noSing) > 1,
                                                   AllREadsPE_Indo_Counts_physeq_noSing)

## Add sequencing depth information after filtering out singletons
SeqDepth_noSingletons <- colSums(otu_table(AllREadsPE_Indo_Counts_physeq_noSing))
sample_data(AllREadsPE_Indo_Counts_physeq_noSing)$SeqDepth_noSingletons <- SeqDepth_noSingletons

## Filter out Viridiplantae 
AllREadsPE_Indo_Counts_physeq_NoPlant <- subset_taxa(AllREadsPE_Indo_Counts_physeq_noSing, (Kingdom!="Viridiplantae"))
AllREadsPE_Indo_Counts_physeq_NoPlant <- prune_taxa(taxa_sums(AllREadsPE_Indo_Counts_physeq_NoPlant) > 0, AllREadsPE_Indo_Counts_physeq_NoPlant)
# add sequencing depth information after filtering out plants
SeqDepth_noViridiplantae = colSums(otu_table(AllREadsPE_Indo_Counts_physeq_NoPlant))
sample_data(AllREadsPE_Indo_Counts_physeq_NoPlant)$SeqDepth_noViridiplantae = SeqDepth_noViridiplantae

## Filter out Chordata
AllREadsPE_Indo_Counts_physeq_noChor <- subset_taxa(AllREadsPE_Indo_Counts_physeq_NoPlant, (Phylum!="Chordata"))
AllREadsPE_Indo_Counts_physeq_noChor <- prune_taxa(taxa_sums(AllREadsPE_Indo_Counts_physeq_noChor) > 0, AllREadsPE_Indo_Counts_physeq_noChor)
# add sequencing depth information after filtering out Metazoa
SeqDepth_noChordata = colSums(otu_table(AllREadsPE_Indo_Counts_physeq_noChor))
sample_data(AllREadsPE_Indo_Counts_physeq_noChor)$SeqDepth_noChordata = SeqDepth_noChordata

## Filter out Metazoa
AllREadsPE_Indo_Counts_physeq_noMeta <- subset_taxa(AllREadsPE_Indo_Counts_physeq_noChor, (Kingdom!="Metazoa"))
AllREadsPE_Indo_Counts_physeq_noMeta <- prune_taxa(taxa_sums(AllREadsPE_Indo_Counts_physeq_noMeta) > 0, AllREadsPE_Indo_Counts_physeq_noMeta)
# add sequencing depth information after filtering out Metazoa
SeqDepth_noMetazoa = colSums(otu_table(AllREadsPE_Indo_Counts_physeq_noMeta))
sample_data(AllREadsPE_Indo_Counts_physeq_noMeta)$SeqDepth_noMetazoa = SeqDepth_noMetazoa

Indo_101BP_PE_postCC <- data.frame(sample_data(AllREadsPE_Indo_Counts_physeq_noMeta))
Indo_101BP_PE_postCC$SampleName <- gsub(".ccm.*","", Indo_101BP_PE_postCC$SampleName)

#umdir <- "/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/QC/ReadsSummary/"
#write.table(Indo_101BP_PE_postCC[,c(1,4:9)], file = paste0(sumdir,"Indo_101BP_PE_postCC.txt"))

# Summarising the data ##########
## Barplot of library sizes
#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/Indo 101BP_PE/libsizes.png"), units="px", width=2500, height=2000, res=400)
ggplot(meta(AllREadsPE_Indo_Counts_physeq_noMeta), aes(SampleName, SeqDepth_noMetazoa)) + geom_bar(stat = "identity", aes(fill = SamplePop)) +
  scale_fill_manual(values = c(KorowaiCol,MentawaiCol,SumbaCol)) + rotate_x_text()
#dev.off()

## Library sizes are highly uneven, with some samples having a  very high library depth (many in the Korowai), and many having a very low library depth.
## Summarize the data and see how many reads lost at each filtering step (on log scale due to large variance in depth between samples). 

FilteringSummary = sample_data(AllREadsPE_Indo_Counts_physeq_noMeta)[,c("SamplePop","SeqDepth_Prefilter","SeqDepth_noSingletons","SeqDepth_noViridiplantae","SeqDepth_noChordata","SeqDepth_noMetazoa")]

FilteringSummary_101BP_PE <- FilteringSummary
rownames(FilteringSummary_101BP_PE) <- sapply(str_split(rownames(FilteringSummary_101BP_PE), ".ccm"),`[`,1)


## melt df and plot
melted_FilteringSummary = melt(FilteringSummary)
ggplot(melted_FilteringSummary, aes(x=variable, y=log(value))) +
  geom_violin(alpha=0.8) + theme_bw() + ylab("Spearman pairwise correlation") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90)) + scale_fill_manual(values = c(KorowaiCol,MentawaiCol,SumbaCol)) +
  geom_boxplot(color="black",width=0.2, alpha = 0.7)

## Most of the reads are removed when removing Chordates.

# Data normalization ##########
## Centered log-ration transformation
## Taxa can be viewed by their relative abundance, however changes in the abundance of one taxon will result in changing the abundance of other taxa.
## One of the ways to handle this is to transform the data using Centered Log Ratio (CLR) transformation. CLR data shows how OTUs behave relative to the per-sample average and is a commonly-used data transformation method in microbiomics.

## CLR transformation is quite sensitive to zeros, so we need to remove as many zeros as we can before performing CLR transformation.
## We'll do this by merging our taxa data at the phylum level, and add an offset of 0.1.
pop_comparison <- AllREadsPE_Indo_Counts_physeq_noMeta %>%
  tax_glom("Phylum")

data_offset <- otu_table(pop_comparison) + 0.1
clr <- clr(data_offset)

## Add in compositions information to new phyloseq object
merged_phylo_counts_clr <- pop_comparison
taxa <- otu_table(clr, taxa_are_rows = TRUE)
otu_table(merged_phylo_counts_clr) <- taxa

#save(merged_phylo_counts_clr, file="/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/AllREadsPE_Indo_Counts_physeq_clr.Rda")

phylum_proportion <- as.data.frame(rowSums(otu_table(pop_comparison)))/sum(rowSums(otu_table(pop_comparison)))

phylum_df <- data.frame(otu_table(pop_comparison))
phylum_md <- data.frame(sample_data(pop_comparison))
phylum_md$SampleName <- str_replace_all(phylum_md$SampleName, "-", ".")

KOR_phyl_df <- phylum_df[,phylum_md[phylum_md$SamplePop == "KOR", "SampleName"]]
KOR_phyl_prop <- data.frame(rowSums(KOR_phyl_df)/sum(rowSums(KOR_phyl_df)))

MTW_phyl_df <- phylum_df[,phylum_md[phylum_md$SamplePop == "MTW", "SampleName"]]
MTW_phyl_prop <- data.frame(rowSums(MTW_phyl_df)/sum(rowSums(MTW_phyl_df)))

SMB_phyl_df <- phylum_df[,phylum_md[phylum_md$SamplePop == "SMB", "SampleName"]]
SMB_phyl_prop <- data.frame(rowSums(SMB_phyl_df)/sum(rowSums(SMB_phyl_df)))


## Also save the data at family level
noMeta_fam <- AllREadsPE_Indo_Counts_physeq_noMeta %>%
  tax_glom("Family")
family_proportion <- as.data.frame(rowSums(otu_table(noMeta_fam)))/sum(rowSums(otu_table(noMeta_fam)))

# Sample clustering ##########
## PCA plots
### Make an ordination plot using euclidean distances
ord <- ordinate(merged_phylo_counts_clr, method = "PCoA", distance = "euclidean")
base_pca <- cbind(ord$vectors, sample_data(merged_phylo_counts_clr))

### Plot by Island
ggplot(base_pca) + geom_point(aes(x=Axis.1, y=Axis.2, color=SamplePop), alpha=1, size=4) + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) + theme(plot.title = element_text(hjust = 0, size = 12))

ggplot(base_pca) + geom_point(aes(x=Axis.3, y=Axis.4, color=SamplePop), alpha=1, size=4) + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) + theme(plot.title = element_text(hjust = 0, size = 12))

ggplot(base_pca) + geom_point(aes(x=Axis.5, y=Axis.6, color=SamplePop), alpha=1, size=4) + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) + theme(plot.title = element_text(hjust = 0, size = 12))

### Sample mainly not by island. We then coor based on pathogenic load:
### Plot by plasmodium load
log10_plas_counts = log10(colSums(otu_table(AllREadsPE_Indo_Counts_physeq_noMeta)[grep("Plasmodiidae",tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Family"])])+1)
sample_data(merged_phylo_counts_clr)[["Plasmodiidae"]] <- log10_plas_counts
ordPlas <- ordinate(merged_phylo_counts_clr, method = "PCoA", distance = "euclidean")
pcs_vase_plas <- cbind(ordPlas$vectors, sample_data(merged_phylo_counts_clr))

Plas_PCs1_2 <- ggplot(pcs_vase_plas) + geom_point(aes(x=Axis.1, y=Axis.2, color=Plasmodiidae, shape=SamplePop), alpha=0.8, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#A03F03"))(10), limits=c(1, 6))
Plas_PCs1_2

ggplot(pcs_vase_plas) + geom_point(aes(x=Axis.3, y=Axis.4, color=Plasmodiidae, shape=SamplePop), alpha=0.8, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#A03F03"))(10), limits=c(1, 6))

ggplot(pcs_vase_plas) + geom_point(aes(x=Axis.5, y=Axis.6, color=Plasmodiidae, shape=SamplePop), alpha=0.8, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#FFCCBB","#A03F03"))(10), limits=c(1, 6))

### Plot by Flavivirus load
log10_flav_counts <- log10(colSums(otu_table(AllREadsPE_Indo_Counts_physeq_noMeta)[grep("Flaviviridae",tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Family"])])+1)

sample_data(merged_phylo_counts_clr)[["Flaviviridae"]] <- log10_flav_counts
ordFlav <- ordinate(merged_phylo_counts_clr, method <- "PCoA", distance = "euclidean")
pcs_vase_flav <- cbind(ordFlav$vectors, sample_data(merged_phylo_counts_clr))

Flav_PCs1_2 <- ggplot(pcs_vase_flav) + geom_point(aes(x=Axis.1, y=Axis.2, color=Flaviviridae, shape=SamplePop), alpha=0.8, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 6))
Flav_PCs1_2

ggplot(pcs_vase_flav) + geom_point(aes(x=Axis.3, y=Axis.4, color=Flaviviridae, shape=SamplePop), alpha=0.8, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 6))

ggplot(pcs_vase_flav) + geom_point(aes(x=Axis.5, y=Axis.6, color=Flaviviridae, shape=SamplePop), alpha=0.8, size=4) + theme(plot.title = element_text(hjust = 0, size = 12)) + scale_colour_gradientn(colours = colorRampPalette(c("#78c679","#006837"))(20), limits=c(1, 6))

### Plot by batch
### We can also see this isn't caused by batch.
plot_ordination(merged_phylo_counts_clr, ord, color="batch", axes = 1:2, label="SampleName") + scale_colour_manual(values=c("#cb54d6","#3d45c4","black")) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
plot_ordination(merged_phylo_counts_clr, ord, color="batch", axes = 2:3, label="SampleName") + scale_colour_manual(values=c("#cb54d6","#3d45c4","black")) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
plot_ordination(merged_phylo_counts_clr, ord, color="batch", axes = 3:4, label="SampleName") + scale_colour_manual(values=c("#cb54d6","#3d45c4","black")) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))
plot_ordination(merged_phylo_counts_clr, ord, color="batch", axes = 4:5, label="SampleName") + scale_colour_manual(values=c("#cb54d6","#3d45c4","black")) + geom_point(aes(), alpha=0.6, size=4) + theme(plot.title = element_text(hjust = 0, size = 12))


## Hierarchical clustering by Euclidean distance
ps_otu <- data.frame(phyloseq::otu_table(merged_phylo_counts_clr))
ps_otu <- t(ps_otu)
bc_dist <- vegan::vegdist(ps_otu, method = "euclidean")
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
### Provide color codes
meta <- data.frame(phyloseq::sample_data(merged_phylo_counts_clr))
colorCode <- c(KOR = KorowaiCol, `MTW` = MentawaiCol, `SMB` = SumbaCol)
labels_colors(ward) <- colorCode[meta$SamplePop][order.dendrogram(ward)]

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/Indo 101BP_PE/phylo_island.png"), units="px", width=2500, height=2000, res=400)
plot(ward, main="Island")
#dev.off()

### Label colours by Plasmodium load
initial <- .bincode(meta$Plasmodiidae, breaks=seq(min(meta$Plasmodiidae, na.rm=T), max(meta$Plasmodiidae, na.rm=T), len = 80),include.lowest = TRUE)
plasmoCol <- colorRampPalette(c("darkblue", "red", "firebrick"))(80)[initial]
labels_colors(ward)=plasmoCol[order.dendrogram(ward)]

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/Indo 101BP_PE/phylo_Plasmodium.png"), units="px", width=2500, height=2000, res=400)
plot(ward, main="Plasmodium")
#dev.off()

### Label colours by Flavivirus load
initial = .bincode(meta$Flaviviridae, breaks=seq(min(meta$Flaviviridae, na.rm=T), max(meta$Flaviviridae, na.rm=T), len = 80),include.lowest = TRUE)
flavoCol <- colorRampPalette(c("darkblue", "#78c679","#006837"))(79)[initial]
labels_colors(ward)=flavoCol[order.dendrogram(ward)]
plot(ward, main="Flaviviridae")

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/Indo 101BP_PE/phylo_Flaviviridae.png"), units="px", width=2500, height=2000, res=400)
plot(ward, main="Flaviviridae")
#dev.off()

# or, to view it in another way
fviz_dend(ward, k = 3,                 # Cut in three groups
          cex = 0.5,                 # label size
          k_colors = c("#616161", "#78c679", "#800026"),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_bw(),     # Change theme
          main = ""
)

### We can see that Flavivirus and Plasmodium have an effect on expression profiles.

# Relative frequency of taxa ##########
## Add a new column containing family names and superkingdom
tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Superkingdom"] <- paste(tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Superkingdom"], tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Family"], sep="_")

tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Superkingdom"] <- gsub("Bacteria_$", "Bacteria_unclassified", tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Superkingdom"])

tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Superkingdom"] <- gsub("Eukaryota_$", "Eukaryota_unclassified", tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Superkingdom"])

tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Superkingdom"] <- gsub("Eukaryota_unclassified", "Eukaryota_unk", tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Superkingdom"])

tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Superkingdom"] <- gsub("Viruses_$", "Viruses_unclassified", tax_table(AllREadsPE_Indo_Counts_physeq_noMeta)[,"Superkingdom"])

## Create circular bar plots
## We'd take top 20, but only 16 phyla left in this dataset after all filtering
aggregated_phyloCounts <- aggregate_top_taxa(AllREadsPE_Indo_Counts_physeq_noMeta, "Superkingdom", top = 20)
## Transform to relative counts
relative_phyloCounts <- microbiome::transform(aggregated_phyloCounts, "compositional")
## Remove weird extra family names added at the end of Superkingdom names
tax_table(relative_phyloCounts)[,"Superkingdom"] <- paste(sapply(strsplit(taxa_names(relative_phyloCounts), "[_.]"), `[`, 1), sapply(strsplit(taxa_names(relative_phyloCounts), "[_.]"), `[`, 2), sep="_")
## Change "Other_NA" to just "Other"
tax_table(relative_phyloCounts)[,"Superkingdom"][grep("Other", taxa_names(relative_phyloCounts))] <- "Other"

data <- as.matrix(as.data.frame(otu_table(relative_phyloCounts)))
data <- as.data.frame(t(data))
data <- data.frame(data[,colnames(data)[1:11]], Eukaryota_unk=data[,"Eukaryota_unk_f"] + data[,"Eukaryota_unk"], data.frame(data[,colnames(data)[14:16]]))
data$group <- as.character(as.matrix(as.data.frame(phyloseq::sample_data(aggregated_phyloCounts)[,"SamplePop"])))
data$individual <- rownames(data)
data$individual <- as.factor(data$individual)
data$group <- as.factor(gsub("MTW", "Mentawai", gsub("SMB", "Sumba", gsub("KOR", "Korowai", data$group))))

tax_rep <- c("Bacteria: Spirochaetota: Borreliaceae", "Bacteria: Actinomycetota: Nocardioidaceae", "Bacteria: Pseudomonadota: Paracoccaceae", "Bacteria: Actinomycetota: Patulibacteraceae", "Bacteria: Pseudomonadota: Sphaerotilaceae", "Bacteria: unclassified", "Eukaryota: Ascomycota: Cladosporiaceae", "Eukaryota: Basidiomycota: Filobasidiaceae", "Eukaryota: Basidiomycota: Malasseziaceae", "Eukaryota: Ascomycota: Ophiocordycipitaceae", "Eukaryota: Apicomplexa: Plasmodiidae", "Eukaryota: unclassified", "Other", "Viruses: Kitrinoviricota: Flaviviridae", "Viruses: unclassified"
)

colnames(data)[1:15] <- tax_rep

## Transform data in a tidy format (long format)
data = data %>% gather(key = "observation", value="value", -c(16,17)) 

## Set a number of 'empty bar' to add at the end of each group
empty_bar=3
nObsType=nlevels(as.factor(data$observation))
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar*nObsType )
data=rbind(data, to_add)
data=data %>% arrange(group, individual)
data$id=rep( seq(1, nrow(data)/nObsType) , each=nObsType)
data$observation <- factor(tax_rep, levels = tax_rep)

# prepare a data frame for base lines
base_data=data %>% 
  group_by(group) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1

# Make the plot
p <- ggplot(data) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.8) +
  # Add a valu=100/75/50/25 lines. I do it at the beginning to make sure barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.25, xend = start, yend = 0.25), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.50, xend = start, yend = 0.50), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.75, xend = start, yend = 0.75), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1, xend = start, yend = 1), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),5), y = c(0, 0.25, 0.50, 0.75, 1), label = c("0", "0.25", "0.50", "0.75", "1") , color="black", size=2 , angle=0, fontface="bold", hjust=1) +
  
  theme_minimal() +
  coord_polar() + 
  geom_text(data=base_data, aes(x = title, y = 1.2, label=group), hjust=c(0,0,0), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE) +
  scale_fill_manual(values=Merged_Palette) +
  theme(
    axis.ticks=element_blank(),
    axis.ticks.length = unit(0, "pt"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.text = element_text(size=6),
    legend.key.size = unit(0.5, 'cm'),
    plot.margin = margin(0,0,0,0,"pt")
  ) + labs(fill = "")

p 

## Plot Figure 1
patho_plt <- (Plas_PCs1_2 | Flav_PCs1_2) + plot_layout(guides = 'collect') & theme(legend.box = "vertical", legend.direction = "vertical", legend.text = element_text(size=6), legend.key.size = unit(0.5, 'cm')) 

#pdf(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/Indo 101BP_PE/nuFig1.pdf"), width=8, height=8)
(p + guide_area() + plot_layout(guides = 'collect', widths = c(2.5,1))) / patho_plt + plot_annotation(tag_levels = 'A') + plot_layout(heights = c(1.3, 1)) & theme(legend.text = element_text(size=8), legend.key.size = unit(0.4, 'cm'), legend.title = element_text(size=9), plot.tag = element_text(face = 'bold'))
#dev.off()

# Differential abundance testing ##########
# We're interested in testing whether the species composition between islands is significantly different.
# Aldex2 corrects for uneven library sizes, so rarefying is not necessary. The only input needed is the data with singletons removed.
# Since Aldex2 works pairwise, the dataset is subsetted into groups of two. 

## Korowai vs Mentawai
KORVsMTW <- subset_samples(pop_comparison, SamplePop != "SMB")
### Remove any 0s from the data
KORVsMTW <- prune_taxa(taxa_sums(KORVsMTW) > 0, KORVsMTW)
KORVsMTW <- prune_taxa(rownames(tax_table(KORVsMTW))[!grepl("unk_", tax_table(KORVsMTW)[,"Phylum"])], KORVsMTW)
taxa_names(KORVsMTW)=make.unique(tax_table(KORVsMTW)[,"Phylum"])

to_remove <- colnames(otu_table(KORVsMTW))[(colSums(otu_table(KORVsMTW))==0)]
KORVsMTW <- prune_samples(!(sample_names(KORVsMTW) %in% to_remove), KORVsMTW)

### Run aldex2
aldex2_KORVsMTW <- ALDEx2::aldex(data.frame(phyloseq::otu_table(KORVsMTW)), phyloseq::sample_data(KORVsMTW)$SamplePop, test="t", effect = TRUE)

sig_aldex2_KORVsMTW <- aldex2_KORVsMTW %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

sig_aldex2_KORVsMTW
aldex2_KORVsMTW <- aldex2_KORVsMTW[order(aldex2_KORVsMTW$wi.eBH),]

KORVsMTW_DA <- ggplot(aldex2_KORVsMTW, aes(x=effect, y=-log10(wi.eBH))) +
  geom_point() + labs(title='Korowai versus Mentawai') +
  xlab("Effect size") + ylab("-log10(P Value)") + scale_color_manual(values=c("#3b59ad","#262424","#b00b0b")) + theme(legend.position = "none", plot.title = element_text(hjust = 1)) + geom_text_repel(data=head(aldex2_KORVsMTW, 3), aes(label=rownames(aldex2_KORVsMTW)[1:3]), size=2.7) + xlim(-1,1) + ylim(0,0.4)

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/Indo 101BP_PE/KORVsMTW_DA.png"), units="px", width=1500, height=1200, res=350)
KORVsMTW_DA
#dev.off()

## Korowai vs Sumba
KORVsSMB=subset_samples(pop_comparison, SamplePop != "MTW")
KORVsSMB <- prune_taxa(taxa_sums(KORVsSMB) > 0, KORVsSMB)
KORVsSMB <- prune_taxa(rownames(tax_table(KORVsSMB))[!grepl("unk_", tax_table(KORVsSMB)[,"Phylum"])], KORVsSMB)
taxa_names(KORVsSMB)=make.unique(tax_table(KORVsSMB)[,"Phylum"])

to_remove <- colnames(otu_table(KORVsSMB))[(colSums(otu_table(KORVsSMB))==0)]
KORVsSMB <- prune_samples(!(sample_names(KORVsSMB) %in% to_remove), KORVsSMB)

aldex2_KORVsSMB <- ALDEx2::aldex(data.frame(phyloseq::otu_table(KORVsSMB)), phyloseq::sample_data(KORVsSMB)$SamplePop, test="t", effect = TRUE)
sig_aldex2_KORVsSMB <- aldex2_KORVsSMB %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

sig_aldex2_KORVsSMB
aldex2_KORVsSMB <- aldex2_KORVsSMB[order(aldex2_KORVsSMB$wi.eBH),]

KORVsSMB_DA <- ggplot(aldex2_KORVsSMB, aes(x=effect, y=-log10(wi.eBH))) +
  geom_point() + labs(title='Korowai versus Sumba') +
  xlab("Effect size") + ylab("-log10(P Value)") + scale_color_manual(values=c("#3b59ad","#262424","#b00b0b")) + theme(legend.position = "none", plot.title = element_text(hjust = 1)) + geom_text_repel(data=head(aldex2_KORVsSMB, 3), aes(label=rownames(aldex2_KORVsSMB)[1:3]), size=2.7) + xlim(-1,1) + ylim(0,0.4)

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/Indo 101BP_PE/KORVsSMB_DA.png"), units="px", width=1500, height=1200, res=350)
KORVsSMB_DA
#dev.off()

## Mentawai vs Sumba
MTWVsSMB=subset_samples(pop_comparison, SamplePop != "KOR")
any(taxa_sums(MTWVsSMB) == 0)
MTWVsSMB <- prune_taxa(taxa_sums(MTWVsSMB) > 0, MTWVsSMB)
MTWVsSMB <- prune_taxa(rownames(tax_table(MTWVsSMB))[!grepl("unk_", tax_table(MTWVsSMB)[,"Phylum"])], MTWVsSMB)
taxa_names(MTWVsSMB)=make.unique(tax_table(MTWVsSMB)[,"Phylum"])

to_remove <- colnames(otu_table(MTWVsSMB))[(colSums(otu_table(MTWVsSMB))==0)]
MTWVsSMB <- prune_samples(!(sample_names(MTWVsSMB) %in% to_remove), MTWVsSMB)

aldex2_MTWVsSMB <- ALDEx2::aldex(data.frame(phyloseq::otu_table(MTWVsSMB)), phyloseq::sample_data(MTWVsSMB)$SamplePop, test="t", effect = TRUE)
sig_aldex2_MTWVsSMB <- aldex2_MTWVsSMB %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

aldex2_MTWVsSMB <- aldex2_MTWVsSMB[order(aldex2_MTWVsSMB$wi.eBH),]

MTWVsSMB_DA <- ggplot(aldex2_MTWVsSMB, aes(x=effect, y=-log(wi.eBH))) +
  geom_point() + labs(title='Mentawai versus Sumba') +
  xlab("Effect size") + ylab("-log10(P Value)") + scale_color_manual(values=c("#3b59ad","#262424","#b00b0b")) + theme(legend.position = "none", plot.title = element_text(hjust = 1)) + geom_text_repel(data=head(aldex2_MTWVsSMB, 5), aes(label=rownames(aldex2_MTWVsSMB)[1:5]), size=2.7) + xlim(-1,1) + ylim(0,0.4)

#png(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/Indo 101BP_PE/MTWVsSMB_DA.png"), units="px", width=1500, height=1200, res=350)
MTWVsSMB_DA
#dev.off()

## There are no statistically significant differentially abundant taxa after BH correction.

# Alpha diversity ##########
# Alpha diversity measures within-sample diversity and looks at how many taxa are observed, as well as how evenly they are distributed. 
# Divnet is a method for estimating within- and between-community diversity in ecosystems where taxa interact via an ecological network. It accounts for differences in sequencing depth and estimates the number of missing species based on the sequence depth and number of rare taxa in the data.
# To use DivNet, you need unsubsampled data without removing singletons.
## Remove Viridiplantae and Metazoa
AllREadsPE_Indo_Counts_physeq_withSingletons <- subset_taxa(AllREadsPE_Indo_Counts_physeq_withSingletons, (Kingdom!="Viridiplantae"))
AllREadsPE_Indo_Counts_physeq_withSingletons <- subset_taxa(AllREadsPE_Indo_Counts_physeq_withSingletons, (Kingdom!="Metazoa"))
## remove any empty rows
AllREadsPE_Indo_Counts_physeq_withSingletons <- prune_taxa(taxa_sums(AllREadsPE_Indo_Counts_physeq_withSingletons) > 0, AllREadsPE_Indo_Counts_physeq_withSingletons)

## DivNet is computationally expensive, and therefore performing at phylum level. DivNet is run without specifying any hypothesis testing.
## Comparing diversity at the phylum level
pop_comparison <- AllREadsPE_Indo_Counts_physeq_withSingletons %>%
  tax_glom("Phylum")
pop_comparison <- prune_taxa(taxa_sums(pop_comparison) > 0, pop_comparison)
## Samples to keep
keep_samples <- names(colSums(otu_table(pop_comparison))[colSums(otu_table(pop_comparison))>1])
pop_comparison <- prune_samples(keep_samples, pop_comparison)

## If we don't change the sample names here from hyphens to periods, we'll get an error later
sample_names(pop_comparison) <- gsub("\\-", ".", sample_names(pop_comparison))

## Run divnet without specifying any hypothesis testing
dv_pop_comparison <- divnet(pop_comparison, ncores = 8, base = pick_base(t(otu_table(pop_comparison)), automatic_cutoff = TRUE))

### Plot the results of Shannon diversity
### A higher Shannon index means higher diversity, whereas a lowed index number means lower diversity.
summary_df_shannon <- as.data.frame(dv_pop_comparison$shannon %>%
                                      summary %>%
                                      add_column("SampleNames" = pop_comparison %>% otu_table %>% sample_names) %>%
                                      add_column("SamplePop" =  pop_comparison %>% sample_data %>% .[,"SamplePop"] %>% as.matrix(.) %>% .[,1] %>% unname(.)))

Indo_Shan <- ggplot(summary_df_shannon, aes(y = estimate, x = SamplePop, fill = SamplePop)) + geom_violin(alpha=0.7) + 
  geom_jitter(height = 0, width = .2) + geom_boxplot(width=0.08, outlier.color = NA) +
  scale_fill_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) + ggtitle("Shannon Diversity") +
  ylab("Estimate of Shannon Diversity") + stat_summary(fun = "median", colour = "darkblue", size = 4,
                                                       geom = "text", aes(label = round(after_stat(y), 3)),
                                                       position = position_nudge(x = -0.35, y=-0.02))

Indo_Shan

plot(dv_pop_comparison$shannon, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol))

### Plot the results of Simpson's diversity
### The higher the value, the lower the diversity. It measures the probability that two individuals randomly selected from a sample will belong to the same species.
### With this index, 0 represents infinite diversity and 1, no diversity.
summary_df_simpson <- as.data.frame(dv_pop_comparison$simpson %>%
                                      summary %>%
                                      add_column("SampleNames" = pop_comparison %>% otu_table %>% sample_names) %>%
                                      add_column("SamplePop" =  pop_comparison %>% sample_data %>% .[,"SamplePop"] %>% as.matrix(.) %>% .[,1] %>% unname(.)))

ggplot(summary_df_simpson, aes(y = estimate, x = SamplePop, fill = SamplePop)) + geom_violin(alpha=0.7) + 
  geom_jitter(height = 0, width = .2) + geom_boxplot(width=0.08, outlier.color = NA) +
  scale_fill_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) + ggtitle("Simpson's Diversity Index") +
  ylab("Estimate of Simpson Diversity")

plot(dv_pop_comparison$simpson, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol))

### Since a larger Simpson index value equates to a lower diversity index, many people find this confusing and not very intuitive.
### Therefore, the inverse Simpson's Index, or 1 - Simpson Index, is also commonly used.
### Subtract the Simpson estimate from one
summary_df_simpson$negEstimate = 1-summary_df_simpson$estimate

### Plot
Indo_InvS <- ggplot(summary_df_simpson, aes(y = negEstimate, x = SamplePop, fill = SamplePop)) + geom_violin(alpha=0.7) + 
  geom_jitter(height = 0, width = .2) + geom_boxplot(width=0.08, outlier.color = NA) +
  scale_fill_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) + ggtitle("Inverse Simpson's Diversity") +
  ylab("Inverse Simpson's Diversity") + stat_summary(fun = "median", colour = "darkblue", size = 4,
                                                     geom = "text", aes(label = round(after_stat(y), 3)),
                                                     position = position_nudge(x = -0.35, y=-0.01))

Indo_InvS
aggregate(summary_df_simpson[, 1], list(summary_df_shannon$SamplePop), median)

#pdf(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/Indo 101BP_PE/nuSFig3_AlphaIndo.pdf"), width=7.5, height=5)
(Indo_Shan | Indo_InvS) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold'), plot.title = element_text(hjust = 1), text = element_text(size=9))
#dev.off()

### Test the hypothesis that the diversity is different between islands
dv_pop_comparison_cov <- pop_comparison %>%
  divnet(X = "SamplePop", ncores = 8, base = pick_base(t(otu_table(pop_comparison)), automatic_cutoff = TRUE))
plot(dv_pop_comparison_cov$shannon, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol))

# Beta diversity ##########
# Beta diversity is a measure of dissimilarity metric between samples to compare differences in species composition.
# It's helpful to know not only how taxonomically/pathogenically rich each sample is, but also to see differences in samples and populations.

## Bray-curtis dissimilarity at the individual sample level
## Bounded between 0 and 1, where 0 means the two sites have the same composition (all species are shared), and 1 means the two sites do not share any species.
bray_est <- simplifyBeta(dv_pop_comparison, pop_comparison, "bray-curtis", "SamplePop")

### Add in group comparisons and plot
bray_est$group=paste(bray_est$Covar1,bray_est$Covar2,sep="_")
islands_BCE <- ggplot(bray_est, aes(x = interaction(Covar1, Covar2), y = beta_est, fill=group)) +
  geom_violin(alpha=0.7) + geom_boxplot(width=0.1) + xlab("Population Comparisons") + 
  theme(legend.position="none") + ggtitle("Bray-Curtis Distance Estimate") +
  ylab("Bray-Curtis Distance") + stat_summary(fun = "mean", colour = "red", size = 4,
                                              geom = "text", aes(label = round(after_stat(y), 3)),
                                              position = position_nudge(x = -0.35, y=-0.05)) + ylim(0,1)

# Plot Figure 2
#pdf(paste0("/Users/muhamad.fachrul/OneDrive - St Vincent's Institute/Projects/Blood microbiome Indo/reScratch/Indo_data/Indo 101BP_PE/nuFigure 2_noBlanks.pdf"), width=7.5, height=5)
(KORVsMTW_DA | KORVsSMB_DA | MTWVsSMB_DA) / islands_BCE + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold'), plot.title = element_text(hjust = 1), text = element_text(size=9))
#dev.off()