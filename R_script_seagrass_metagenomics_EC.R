### Scripts for working with SQMtools and other R packages to analyze metagenomics data ###
## July 2026
## by Evan Cronmiller

# loading required packages
library(dplyr)
library(ggplot2)
library(tibble)
library(SQMtools)
library(ggpattern)
library(microeco)
library(indicspecies)
library(viridis)
library(MicrobiomeStat)
library(WGCNA)
library(ggtree)
library(agricolae)
library(ggpubr)
library(vegan)
library(RColorBrewer)
library(FSA)
library(rcompanion)
library(magick)
library(SBGNview)
library(DESeq2)

## quick look at some Squeezemeta output stats
sqmstats <- read.delim ("C:/Users/CRONMILLERE/Documents/research_related/COOGER Metagenomics/read_mapping_stats.txt")
# total reads
mean(sqmstats$Total.reads)
range(sqmstats$Total.reads)
# reads mapped
mean(sqmstats$Mapping.perc)
range(sqmstats$Mapping.perc)

###### for SQMtools analysis of data from SqueezeMeta pipeline ######
# loading in data from SqueezeMeta (zip file)
sgsqm = loadSQM("D:/Seagrass_Meta/Seagrass.zip",
                      load_sequences = FALSE)  #not loading seqs right now because too much memory

# bringing in metadata file for seagrass project samples
sgmetadat <- read.csv("field_meta_data.csv",
                      blank.lines.skip = TRUE)

# adding to the metadata the sample groups identified via 16S analysis
sg16sgroups <- as.data.frame(read.table("sample_groups_from_16s.tsv",
                                        header = TRUE,
                                        col.names = c("ID", "group_16s"))) # changing the column names to enable dataframe matching
# modifying sample IDs to match sample IDs in metadata (removing S prefix) and making numeric
sg16sgroups$ID <- as.numeric(gsub("^S", "", sg16sgroups$ID)) 

# combining dataframes by sample ID to create single metadata file to be used for metagenome analyses
# this will exclude samples from the metadata that were not used for metagenomic sequencing
sgmetadatcomb <- inner_join(sgmetadat, sg16sgroups, by = "ID")

# bringing in environmental data from E. Cormier to add to our metadata
sgenvdat <- readxl::read_xlsx("OPP-ERS-ES-2025-FullEnviroData.xlsx", sheet = "OPP-ERS-ES-2025-FullEnviroData")
sgmetadatcomb <- merge(sgmetadatcomb, sgenvdat, by = c("Site", "Quadrat"))
# quickly looking at combined metadata to see if environmental data are appropriately combined
ggplot(sgmetadatcomb, aes(x = Site)) + 
  geom_point(aes(y = n.veg.shoots.per.m2, color = "shoots_m2")) +
  geom_point(aes(y = Per.OM, color = "organic_matter")) +
  geom_point(aes(y = Per.mud.silt, color = "percent_mud"))


## checking BINs and contigs for data from SqueezeMeta pipeline
# looking at BIN statistics from SqueezeMeta
binstats <- data.frame(row.names = row.names(sgsqm$bins$table), sgsqm$bins$table$Completeness, sgsqm$bins$table$Contamination) # subsetting the relevant BIN statistics
binstats <- rownames_to_column(binstats, var = "bin_name")   # preparing data for plotting
binstats$bin_name <- factor(binstats$bin_name, levels = unique(binstats$bin_name))
# plotting bin statistics for completeness and contamination  
ggplot(binstats, aes(x= bin_name)) +
  geom_point(aes(y = sgsqm.bins.table.Completeness, color = "Completeness")) +
  geom_point(aes(y = sgsqm.bins.table.Contamination, color = "Contamination")) +
  theme(axis.text.x = element_blank(), legend.title = element_blank()) +
  ylab("Percent") +
  xlab("BIN")

# how many BINs remain if we apply cutoffs for completeness and contamination?
binstatsfilt <- subset(binstats, sgsqm.bins.table.Completeness > 70 & sgsqm.bins.table.Contamination < 10)
print(nrow(binstatsfilt))

### visualizing metagenomics data using SQMtools ###
## barplot of taxonomy adbundances by taxonomic rank 
plotBars(sgsqm$taxa$superkingdom$abund)  # for superkingdom

plotBars(sgsqm$taxa$phylum$abund)  # for phylum (too many phyla to plot)

# creating a list object with the sample groups to enable this vizualization with samples grouped together
groupsite <- list(PH = c(518514, 518515, 518517), FishPed = c(518524, 518526, 518528), PJ = c(518535, 518540, 518542), CresB = c(518545, 518547, 518549), L3F = c(518555, 518557, 518559), CblSH = c(518563, 518565, 518567), CblD = c(518568, 518570, 518571), Crch = c(518585, 518587, 518591), SambD = c(518595, 518597, 518601), SambSH = c(518605, 518607, 518609), SacrD = c(518615, 518617, 518619), SacrSH = c(518623, 518627, 518629), TH = c(518633, 518637, 518641), PadH = c(518645, 518647, 518651))
group16s <- list(groupA = c(518514, 518515, 518517,518524, 518526, 518528,518535, 518540, 518542), groupB = c(518545, 518547, 518549,518555, 518557, 518559), groupC = c(518595, 518597, 518601, 518605, 518607, 518609, 518633, 518637, 518641, 518645, 518647, 518651), groupD = c(518563, 518565, 518567,518568, 518570, 518571,518585, 518587, 518591,518615, 518617, 518619, 518623, 518627, 518629))     # groupings from 16S data analysis (Alice)

# trying again barplots, playing around with combining by groups and removing unclassified/unmapped hits
plotBars(sgsqm$taxa$superkingdom$abund, metadata_groups = groupsite)
plotBars(sgsqm$taxa$superkingdom$percent, metadata_groups = groupsite,
         label_y = "Percent",
         label_x = "Samples (grouped by site)")

plotBars(sgsqm$taxa$superkingdom$percent, metadata_groups = group16s,  # trying again with 16s groupings
         label_y = "Percent",
         label_x = "Samples (grouped by 16s data)")


## plotting the most abundant taxa
# at phyla level

plotTaxonomy(sgsqm,   # data source
             rank = "phylum",   # taxonomic level to pull data from
             count = "abund",   # type of metric/data to plot (i.e. abundance vs. percentage)
             N = 15)   #number of top hits to include

plotTaxonomy(sgsqm,   # trying with class data 
             rank = "class",   # taxonomic level to pull data from
             count = "abund",   # type of metric/data to plot (i.e. abundance vs. percentage)
             N = 15)   #number of top hits to include

# trying again but with additional parameters, by percentage (easier to compare visually) and without unmapped reads
plotTaxonomy(sgsqm,  
             rank = "class", 
             count = "percent",   
             N = 15,
             metadata_groups = groupsite,   # grouped by site
             ignore_unmapped = TRUE,   # removing unmapped reads
             ignore_unclassified = TRUE,
             nocds = "treat_as_unclassified")  # removing unclassified reads

plotTaxonomy(sgsqm,  
             rank = "class", 
             count = "percent",   
             N = 15,
             metadata_groups = group16s,   # grouped by 16s result
             ignore_unmapped = TRUE,   # removing unmapped reads
             ignore_unclassified = TRUE,
             nocds = "treat_as_unclassified")  # removing unclassified reads

plotTaxonomy(sgsqm,  
             rank = "genus", 
             count = "percent",   
             N = 15,
             metadata_groups = groupsite,   # grouped by 16s result
             ignore_unmapped = TRUE,   # removing unmapped reads
             ignore_unclassified = TRUE,
             nocds = "treat_as_unclassified")  # removing unclassified reads


# trying again but with changed color palette
colpal <- magma(20)  # viridis colour palette
plotTaxonomy(sgsqm,  
             rank = "class", 
             count = "percent",   
             N = 20,
             metadata_groups = groupsite,   # grouped by site
             ignore_unmapped = TRUE,   # removing unmapped reads
             ignore_unclassified = TRUE,
             nocds = "treat_as_unclassified",
             color = colpal)    # hard to see difference... need a randomized colour palette


## plotting functional data (i.e. KEGG, COG, PFAM)
plotFunctions(sgsqm,
              fun_level = "KEGG",   # choosing KEGG data
              count = "tpm",  # using tpm datasets
              N = 25,         # top 25 functions
              metadata_groups = groupsite)  # number of most abundant functional units to plot

plotFunctions(sgsqm,
              fun_level = "KEGG",   # choosing KEGG data
              count = "tpm",  # using tpm datasets
              N = 25,         # top 25 functions
              metadata_groups = group16s)  # grouped by 16s data

# with cpm as hit counting metric
plotFunctions(sgsqm,
              fun_level = "KEGG",   # choosing KEGG data
              count = "cpm",  # using cpm datasets
              N = 25,         # top 25 functions
              metadata_groups = groupsite)  # number of most abundant functional units to plot

# trying again but with COG data
plotFunctions(sgsqm,
              fun_level = "COG",   # choosing COG data
              count = "tpm",  # using tpm datasets
              N = 25,
              gradient_col = c("ghostwhite", "darkgreen"),
              metadata_groups = groupsite)  # number of most abundant functional units to plot

plotFunctions(sgsqm,
              fun_level = "COG",   # choosing COG data
              count = "tpm",  # using tpm datasets
              N = 25,
              gradient_col = c("ghostwhite", "darkgreen"),
              metadata_groups = group16s)  # grouped by 16s data


# trying again but with PFAM data
plotFunctions(sgsqm,
              fun_level = "PFAM",   # choosing PFAM data
              count = "cpm",  # using cpm datasets
              N = 25,         # top 25 functions
              metadata_groups = groupsite)  # number of most abundant functional units to plot


### trying pathway function in SQMtools which is a wrapper for the R package "pathway" 
# trying for KEGG pathway K00920, which is sulfate metabolism
# testing with fold change difference in median tpm across samples within the two groups  

# making groups for pathway testing
fishped <- c("518524", "518526", "518528")
sambd <- c("518595", "518597", "518601")

# trying for KEGG pathway K00920, which is sulfur metabolism, comparing two groups with expected differences
exportPathway(sgsqm,
              pathway_id = "00920",
              count = "tpm",
              fold_change_groups = list(fishped, sambd),
              output_suffix = "sulf_path",
              output_dir = "C:/Users/CRONMILLERE/Documents/research_related/COOGER Metagenomics/Figures")
# trying for KEGG pathway K00680, which is methane metabolism
exportPathway(sgsqm,
              pathway_id = "00680",
              count = "tpm",
              fold_change_groups = list(fishped, sambd),
              output_suffix = "methane_meta__fished_samb_FoldChange",
              output_dir = "C:/Users/CRONMILLERE/Documents/research_related/COOGER Metagenomics/Figures")

# trying once more for pathway K00624, polycyclice aromatic hydrocarbon degradation
exportPathway(sgsqm,
              pathway_id = "00624",
              count = "tpm",
              fold_change_groups = list(fishped, sambd),
              output_suffix = "pahs_meta__fished_samb_FoldChange",
              output_dir = "C:/Users/CRONMILLERE/Documents/research_related/COOGER Metagenomics/Figures")


####### trying export data to be used with microeco package ##########
# making metadata dataframe have rows labelled by sample ID to match format required by microeco
rownames(sgmetadatcomb) <- sgmetadatcomb$ID

# trying a microtable for BINs with abundance as measurement
micbinabund <- SQM_to_microeco(sgsqm,
                               features = "bins",
                               count = "abund",
                               md = sgmetadatcomb)

# trying to create another microtable for feature KEGG
micfuncttpm <- SQM_to_microeco(sgsqm,
                                   features = "KEGG",
                                   count = "tpm",
                                   md = sgmetadatcomb)

# trying basic bar plot to see if conversion to microtable format worked
abundnew <- trans_abund$new(dataset = micbinabund,
                      taxrank = "Class",
                      ntaxa = 10)
abundnew$plot_bar(others_color = "grey70",
            xtext_keep = F, 
            legend_text_italic = FALSE,
            clustering_plot = TRUE)

# creating another microtable with BIN data, measured by cpm, and with unmapped  (doesn't matter anyways)
micbincpm <- SQM_to_microeco(sgsqm,
                               features = "bins",
                               count = "cpm",
                               md = sgmetadatcomb,
                             no_partial_classifications = F,
                             ignore_unclassified = F,
                             ignore_unmapped = F)

# creating another microtable with BIN data, measured by abundance
micbinabund <- SQM_to_microeco(sgsqm,
                               features = "bins",
                               count = "abund",
                               md = sgmetadatcomb,
                               no_partial_classifications = F,
                               ignore_unclassified = F,
                               ignore_unmapped = F)

# creating another microtable with BIN data, measured by abundance, and no unclassified/unmapped hits
micbinabundnounmap <- SQM_to_microeco(sgsqm,
                                     features = "bins",
                                     count = "abund",
                                     md = sgmetadatcomb,
                                     no_partial_classifications = F,
                                     ignore_unclassified = F,
                                     ignore_unmapped = T)

# creating contig table, although much of the useful info here can be taken from the taxa tables...
miccontigsabundnopart <- SQM_to_microeco(sgsqm,
                                     features = "contigs",
                                     count = "abund",
                                     md = sgmetadatcomb,
                                     no_partial_classifications = F,
                                     ignore_unclassified = F,
                                     ignore_unmapped = F)

# creating taxa microtable (i.e. contig hits)
mictaxaabundnounmap <- SQM_to_microeco(sgsqm,
                                         features = "species",
                                         count = "abund",
                                         md = sgmetadatcomb,
                                         no_partial_classifications = FALSE,
                                         ignore_unclassified = TRUE,
                                         ignore_unmapped = TRUE)
tidy_taxonomy(mictaxaabundnounmap, add_prefix = TRUE) #cleaning up taxonomic labels for visualizations

# and again with taxa table but including unmapped hits
mictaxaabund <- SQM_to_microeco(sgsqm,
                                       features = "species",
                                       count = "abund",
                                       md = sgmetadatcomb,
                                       no_partial_classifications = FALSE,
                                       ignore_unclassified = TRUE,
                                       ignore_unmapped = FALSE)
tidy_taxonomy(mictaxaabundnounmap, add_prefix = FALSE) #cleaning up taxonomic labels for visualizations

# trying a heatmap by sample to look at taxa abundances, first from "BINs" SQM table
transbingenus <- trans_abund$new(dataset = micbinabund, taxrank = "Genus", ntaxa = 40)  # choosing top 40 abundant genera
transbingenushmap <- transbingenus$plot_heatmap(facet = "Site", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
transbingenushmap
# and again BIN at phylum level
transbinphylum <- trans_abund$new(dataset = micbinabund, taxrank = "Phylum", ntaxa = 40)  # c
transbinphylumhmap <- transbinphylum$plot_heatmap(facet = "Site", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
transbinphylumhmap
# looking at taxa abundance from SQM taxa tables (i.e. hits to contigs)
# first looking at Phylum level
transtaxphylum <- trans_abund$new(dataset = mictaxaabundnounmap, taxrank = "Phylum", ntaxa = 40)
transtaxphylumhmap <- transtaxphylum$plot_heatmap(facet = "Site", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
transtaxphylumhmap   # data also includes non-prokaryotes!
# now at Genus level
transtaxgenus <- trans_abund$new(dataset = mictaxaabundnounmap, taxrank = "Genus", ntaxa = 40)
transtaxgenushmap <- transtaxgenus$plot_heatmap(facet = "Site", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
transtaxgenushmap   # really not good enough resolution at gennus level, makes viz look bad
# best resolution is at the Class level
transtaxclass <- trans_abund$new(dataset = mictaxaabundnounmap, taxrank = "Class", ntaxa = 40)
transtaxclasshmap <- transtaxclass$plot_heatmap(facet = "Site", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
transtaxclasshmap   # really not good enough resolution at gennus level, makes viz look bad

# prepare for plotting of cpm data 
cpmnew <- trans_abund$new(dataset = micbincpm,
                            taxrank = "Class",
                            ntaxa = 10)
cpmnew$plot_bar(others_color = "grey70",
                  xtext_keep = T, 
                  legend_text_italic = FALSE,
                  clustering_plot = TRUE)  #including dendrogram of sample clustering

# trying again with abundance data to compare with cpm
abundnewnopart <- trans_abund$new(dataset = micbinabund,
                                taxrank = "Class",
                                ntaxa = 10)
abundnewnopart$plot_bar(others_color = "grey70",
                      xtext_keep = T, 
                      legend_text_italic = FALSE,
                      clustering_plot = TRUE)   
# result is very similar, but includes no BIN reads, while cpm measure excludes this

# preparing abundance data averaged by sampling site
cpmnewgroup <- trans_abund$new(dataset = micbincpm,
                                  taxrank = "Class",
                                  ntaxa = 10,
                                  groupmean = "Site")
# plotting
cpmnewgroup$plot_bar(others_color = "grey70",
                        xtext_keep = T, 
                        legend_text_italic = F)


#### trying diversity statistics using microeco ####

# making 14 value colour palette for plots with not enough default colours
dark2 <- brewer.pal(8, "Dark2")
set1 <- brewer.pal(9, "Set1")[c(1:5,7)]  
extended_palette <- c(dark2, set1)
extended_palette <- extended_palette[1:14]

## calculating alpha diversity of datasets
# testing alpha diversity
# first trying with BIN data, including unmapped hits
binabundalpha <- trans_alpha$new(dataset = micbinabund, group = "Site")  # need to use abundances for alpha div. as cpm is not appropriate
head(binabundalpha$data_stat)  #Chao1 is the exact same for all samples... BIN data collapses differences by removing certain taxa...
# trying some comparative stats between sampling sites for species composition with BINs
# first Kruskal-Wallis
binabundalpha$cal_diff(method = "KW")   # nothing is significant... again BIN data seems to remove a lot of naunce
# trying ANOVA
binabundalpha$cal_diff(method = "anova")  # some sampling sites do group together... but not necessarily as expected
# plotting Chao1
binabundalpha$plot_alpha(measure = "Chao1")   # this doesn't work because the Chao1 is the same for all samples with BIN data
# trying again with a different alpha statistic
binabundalpha$plot_alpha(measure = "Shannon", add = "jitter") 
#*** kruskal-Wallis test is probably most appropriate??


# trying again with BIN data, but grouped by 16s groupings
binabundalpha16s <- trans_alpha$new(dataset = micbinabund, group = "group_16s")  # need to use abundances for alpha div. as cpm is not appropriate
head(binabundalpha16s$data_stat)  #Chao1 is the exact same for all samples... BIN data collapses differences by removing certain taxa...
# trying some comparative stats between sampling sites for species composition with BINs
# first Kruskal-Wallis
binabundalpha16s$cal_diff(method = "KW")   
# trying ANOVA
binabundalpha16s$cal_diff(method = "anova")  
# plotting Chao1
binabundalpha16s$plot_alpha(measure = "Chao1")   # this doesn't work because the Chao1 is the same for all samples with BIN data
# trying again with a different alpha statistic
binabundalpha16s$plot_alpha(measure = "Shannon", add = "jitter") 

# now trying with BIN data, excluding unmapped hits
binabundnounmapalpha <- trans_alpha$new(dataset = micbinabundnounmap, group = "Site")  
head(binabundnounmapalpha$data_stat)  
binabundnounmapalpha$cal_diff(method = "KW")   
binabundnounmapalpha$plot_alpha(measure = "Shannon", add = "jitter") 

# trying again with BIN data, no unmapped BIN hits included, but grouped by 16s groupings
binabundnounmapalpha16s <- trans_alpha$new(dataset = micbinabundnounmap, group = "group_16s")  # need to use abundances for alpha div. as cpm is not appropriate
head(binabundnounmapalpha16s$data_stat)  #Chao1 is the exact same for all samples... BIN data collapses differences by removing certain taxa...
binabundnounmapalpha16s$cal_diff(method = "KW")   
binabundnounmapalpha16s$plot_alpha(measure = "Shannon", add = "jitter") 


# now trying instead with contig-derived taxa data  (***note taxa data includes metazoans***)
taxaalpha <- trans_alpha$new(dataset = mictaxaabund, group = "Site")
head(taxaalpha$data_stat)   # better variability between sites and samples with the taxa expanded dataset compared to BINs
taxaalpha$cal_diff(method = "KW_dunn") 
taxaalpha$plot_alpha(measure = "Shannon", add = "jitter", xtext_size = 10)

# trying again with contig-derived taxa data, but grouped by 16s results
taxaalpha16s <- trans_alpha$new(dataset = mictaxaabund, group = "group_16s")
head(taxaalpha16s$data_stat)
taxaalpha16s$cal_diff(method = "KW_dunn") 
taxaalpha16s$plot_alpha(measure = "Shannon", add = "jitter")

# now trying instead with contig-derived taxa data, excluding unmapped reads (***note taxa data includes metazoans***)
taxanounmapalpha <- trans_alpha$new(dataset = mictaxaabundnounmap, group = "Site")
head(taxanounmapalpha$data_stat)   # better variability between sites and samples with the taxa expanded dataset compared to BINs
taxanounmapalpha$cal_diff(method = "KW_dunn") 
taxanounmapalpha$plot_alpha(measure = "Shannon", add = "jitter", xtext_size = 10)

# trying again with contig-derived taxa data, excluding unmapped reads, but grouped by 16s results
taxanounmapalpha16s <- trans_alpha$new(dataset = mictaxaabundnounmap, group = "group_16s")
head(taxanounmapalpha16s$data_stat)
taxanounmapalpha16s$cal_diff(method = "KW_dunn") 
taxanounmapalpha16s$plot_alpha(measure = "Shannon", add = "jitter")


# trying again with differential test of groups (site and 16s)
taxadiffalpha <- trans_alpha$new(dataset = mictaxaabundnounmap, group = "Site", by_group = "group_16s")
taxadiffalpha$cal_diff(method = "KW")
taxadiffalpha$plot_alpha(measure = "Shannon")


#######

## testing beta diversity
# first with BIN abundance data
micbinabund$cal_betadiv()
# creating trans_beta object with Bray-Curtis dissimilarity measure as calculated by cal_betadiv
transbetabinabund <- trans_beta$new(dataset = micbinabund, group = "Site", measure = "bray")
# calculating the ordination with PCoA metric to enable plotting of bray dissimilarity
transbetabinabund$cal_ordination(method = "PCA")  # PCA is more suitable for abundance metrics (as opposed to cpm)
# plotting the ordination, with point colour and symbol by sampling site
transbetabinabund$plot_ordination(plot_color = "Site", plot_shape = "Site",plot_type = "point", point_size = 4, point_alpha = 1, color_values = extended_palette)

# trying again with BIN cpm data  (maybe better because is a normalized value and excludes "no BIN" read data)
micbincpm$cal_betadiv()
transbetabincpm <- trans_beta$new(dataset = micbincpm, group = "Site", measure = "bray") # Jaccard dissimilarity doesn't work with BIN data 
transbetabincpm$cal_ordination(method = "PCoA")
transbetabincpm$plot_ordination(plot_color = "Site", plot_shape = "Site",plot_type = c("point"), point_size = 4, point_alpha = 1, color_values = extended_palette)
# now trying same dataset and bray-curtis dissimilarity, but with NMDS ordination
transbetabincpm$cal_ordination(method = "NMDS")
transbetabincpm$plot_ordination(plot_color = "Site", plot_shape = "Site",plot_type = c("point"), point_size = 4, point_alpha = 1, color_values = extended_palette)

# looking at this same data, just organized by 16s groups
transbetabincpm16s <- trans_beta$new(dataset = micbincpm, group = "group_16s", measure = "bray")
transbetabincpm16s$cal_ordination(method = "PCoA")
transbetabincpm16s$plot_ordination(plot_color = "group_16s", plot_shape = "group_16s",plot_type = c("point", "chull"), point_size = 4, point_alpha = 1)
# and with NMDS
transbetabincpm16s$cal_ordination(method = "NMDS")
transbetabincpm16s$plot_ordination(plot_color = "group_16s", plot_shape = "group_16s",plot_type = c("point", "ellipse"), point_size = 4, point_alpha = 1)


# now with taxa abundance data (excluding unmapped reads)
mictaxaabundnounmap$cal_betadiv()
transbetataxaabund <- trans_beta$new(dataset = mictaxaabundnounmap, group = "Site", measure = "bray")
transbetataxaabund$cal_ordination(method = "PCA")  # PCA or DCA is more suitable for abundance metrics (as opposed to cpm)
transbetataxaabund$plot_ordination(plot_color = "Site", plot_shape = "Site",plot_type = "point", point_size = 4, point_alpha = 1, color_values = extended_palette)
# again, but grouped by 16s group
transbetataxaabund16s <- trans_beta$new(dataset = mictaxaabundnounmap, group = "group_16s", measure = "bray")
transbetataxaabund16s$cal_ordination(method = "PCA")  # PCA or DCA is more suitable for abundance metrics (as opposed to cpm)
transbetataxaabund16s$plot_ordination(plot_color = "group_16s", plot_shape = "group_16s",plot_type = "point", point_size = 4, point_alpha = 1, color_values = extended_palette)


## Trying Jaccard dissimilarity matrix with taxa abundance data
transbetataxaabundjaccard <- trans_beta$new(dataset = mictaxaabundnounmap, group = "Site", measure = "jaccard")
transbetataxaabundjaccard$cal_ordination(method = "NMDS")  # PCA or DCA is more suitable for abundance metrics (as opposed to cpm)
transbetataxaabundjaccard$plot_ordination(plot_color = "Site", plot_shape = "Site",plot_type = "point", point_size = 4, point_alpha = 1, color_values = extended_palette)
 ## looks very similar to Bray-Curtis results! Don't need to do both, I think, just Bray-Curtis should be fine


### exploring trans_diff class from microeco to look at differential abundance ###

## Running trans_diff with taxa data (contigs) to see which taxa are differentially abundant at between sites
# requires rownames to be sample IDs, so changing that here
rownames(mictaxaabundnounmap$sample_table) <- mictaxaabundnounmap$sample_table$ID
## creating trans_diff object using LEfSe model for differential abundance (using taxa abundance data)
testdiff <- trans_diff$new(dataset = mictaxaabundnounmap, method = "lefse", group = "Site", alpha = 0.01)

# plotting trans diff LEfSe results
testdiff$plot_diff_bar(threshold = 4.0, add_sig = T)
testdiff$plot_diff_abund(use_number = 1:20, coord_flip = F, simplify_names = T)
testdiff$plot_diff_cladogram()

# trying again with class taxnomic rank as cutoff
taxaabundnounmapclassdiff <- trans_diff$new(dataset = mictaxaabundnounmap, method = "lefse", group = "Site", alpha = 0.01, taxa_level = "Class")
taxaabundnounmapclassdiff$plot_diff_bar(use_number = 2:51, add_sig = F)
taxaabundnounmapclassdiff$plot_diff_cladogram()

# trying again, but with tidied taxonomy and grouped by 16s result
taxaabundnounmap16sdiff <- trans_diff$new(dataset = mictaxaabundnounmap, method = "lefse", group = "group_16s", alpha = 0.01)
taxaabundnounmap16sdiff$plot_diff_bar(threshold = 3.5, add_sig = T, keep_prefix = T)
taxaabundnounmap16sdiff$plot_diff_abund()
taxaabundnounmap16sdiff$plot_diff_cladogram()

# again with the LEfSe model for BIN data, using abundances
binabundnounmapdiff <- trans_diff$new(dataset = micbinabundnounmap, method = "lefse", group = "Site", alpha = 0.01)
binabundnounmapdiff$plot_diff_bar(threshold = 3.5, add_sig = T)

binabundnounmap16sdiff <- trans_diff$new(dataset = micbinabundnounmap, method = "lefse", group = "group_16s", alpha = 0.01)
binabundnounmap16sdiff$plot_diff_bar(threshold = 3.2, add_sig = T)
# now again with BIN data, but cpm
bincpmnounmapdiff <- trans_diff$new(dataset = micbincpm, method = "lefse", group = "Site", alpha = 0.01)
bincpmnounmapdiff$plot_diff_bar(threshold = 3.8, add_sig = T)
# cpm data closely agrees with abundance, but not exactly


testbincpmdiff <- trans_diff$new(dataset = micbincpm, method = "lefse", group = "Site", alpha = 0.01, taxa_level = "Family")
testbincpmdiff$plot_diff_bar(threshold = 3, add_sig = T)


# testing cladogram with taxa data (which contains the taxonomic rank labels required by cladogram)
t1 <- trans_diff$new(dataset = mictaxaabundnounmap, method = "lefse", group = "group_16s", alpha = 0.01)
t1$res_diff[1:5,] # checking
t1$plot_diff_cladogram(use_taxa_num = 100, use_feature_num = 30, clade_label_level = 5, group_order = c("groupA", "groupB", "groupC", "groupD"), clade_label_size = 1)
# the labels don't show up on the cladogram... not sure how to fix this as its a problem with the formulas/code


## now trying other trans_diff models for differential analysis

binabundnounmapdeseqdiff <- trans_diff$new(dataset = micbinabundnounmap, method = "DESeq2", group = "Site", alpha = 0.01)

binabundnounmapdeseqphylumdiff <- trans_diff$new(dataset = micbinabundnounmap, method = "DESeq2", group = "Site", alpha = 0.01, taxa_level = "Phylum")
binabundnounmapdeseqphylumdiff$plot_diff_bar(heatmap_cell = "log2FoldChange", heatmap_sig = "Significance", heatmap_lab_fill = "log2FoldChange")

binabundnounmapkwdiff <- trans_diff$new(dataset = micbinabundnounmap, method = "KW", group = "Site")

binabundnounmaplindadiff <- trans_diff$new(dataset = micbinabundnounmap, method = "linda", group = "Site", filter_thres = 0.001)


##### now exploring functional diversity with microeco #####

## trying FAPROTAX analysis 
# to work FAPROTAX with trans_func, need to change domain to kingdom (even though not technically correct)
micbinabundnounmapforfun = micbinabundnounmap # cloning the bin abundance dataset
micbinabundnounmapforfun$tax_table$Kingdom <- micbinabundnounmapforfun$tax_table$Domain # copying the domain taxa column to "kingdom"

t2 <- trans_func$new(dataset = micbinabundnounmapforfun)
t2$cal_func(prok_database = "FAPROTAX") # applying FAPROTAX analysis
t2$res_func[1:5, 1:3]
t2$cal_func_FR(abundance_weighted = FALSE, perc = FALSE, adj_tax = FALSE)
t2$res_func_FR[1:5, 1:3]

mictaxaabundnounmapforfun = mictaxaabundnounmap
mictaxaabundnounmapforfun$tax_table$Kingdom <- mictaxaabundnounmapforfun$tax_table$Domain
taxaabundfunc <- trans_func$new(dataset = mictaxaabundnounmapforfun)
taxaabundfunc$cal_func(prok_database = "FAPROTAX")
taxaabundfunc$res_func[1:20, 1:5]
taxaabundfunc$cal_func_FR(abundance_weighted = TRUE, perc = TRUE, adj_tax = FALSE)
taxaabundfunc$res_func_FR[1:20, 1:5]
taxaabundfunc$plot_func_FR()

testdatforfun = mictaxaabund
testdatforfun$tax_table$Kingdom <- testdatforfun$tax_table$Domain
testdatfunc <- trans_func$new(dataset = testdatforfun)
testdatfunc$cal_func(prok_database = "FAPROTAX")
testdatfunc$res_func[1:20, 1:5]
testdatfunc$cal_func_FR(abundance_weighted = TRUE, perc = TRUE, adj_tax = FALSE)
testdatfunc$res_func_FR[1:20, 1:5]
testdatfunc$plot_func_FR()


taxaabundfunc$

tmp_mt <- clone(mictaxaabundnounmap)

taxaabundfuncotu <- as.data.frame(t(taxaabundfunc$res_func_FR), check.names = FALSE)

tmp_mt$taxa_abund$func <- taxaabundfuncotu

t4 <- trans_diff$new(dataset = tmp_mt, method = "anova", group = "Site", taxa_level = "func")

t4$plot_diff_abund(add_sig = T) + ggplot2::ylab("Relative abundance (%)")


####*** need to try faprotax with renaming superkingdom to kingdom



## looking at ORF data 
micorftpm <- SQM_to_microeco(sgsqm,
                             features = "orfs",
                             count = "tpm",
                             md = sgmetadatcomb)

micorftpm2 <- SQM_to_microeco(sgsqm,
                             features = "orfs",
                             count = "tpm",
                             ignore_unclassified = T,
                             md = sgmetadatcomb)








###############################
#### indicspecies analysis ####
###############################

## looking first at BIN abundances
# creating the matrix/dataframe format required by multipatt function, taken from microtable
 

indval = multipatt(binabundindic, sitegroups)


## looking now at contig/taxa abundances
# creating the matrix/dataframe format required by multipatt function, taken from microtable
taxaabundindic <- t(mictaxaabundnounmap$otu_table)
# creating the group vector required (site)
sitegroups <- as.vector(sgmetadatcomb$Site) # taking values from metadata
names(sitegroups) <- rownames(sgmetadatcomb)
# running multipatt
indictaxa <- multipatt(taxaabundindic, sitegroups, control = how(nperm=100)) ###***need to run this again with HPC w/ 999 permutations
summary(indictaxa)
summary(indictaxa, alpha=0.01)
summary(indictaxa, indvalcomp=TRUE)

# trying this again but with no multi-site group comparisons
indictaxanomulti <- multipatt(taxaabundindic, sitegroups, control = how(nperm=999), duleg = TRUE)
summary(indictaxanomulti, alpha = 0.01)


summary(indictaxanomulti, indvalcom=TRUE)
