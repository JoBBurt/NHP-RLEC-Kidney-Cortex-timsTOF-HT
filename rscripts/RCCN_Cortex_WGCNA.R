
library(WGCNA)

setwd("~/Desktop/RCCN_Kidney/Cortex/Batch_1and2/")

metabolite_counts <- my_pro


counts <- metabolite_counts

#rownames(counts) <- counts$PG.UniProtIds

#counts$PG.UniProtIds <- NULL

counts <- t(counts)

gsg = goodSamplesGenes(counts, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(counts)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(counts)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  counts = counts[gsg$goodSamples, gsg$goodGenes]
}





sampleTree = hclust(dist(counts), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


nGenes = ncol(counts)
nSamples = nrow(counts)



# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=50, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(counts, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



allowWGCNAThreads(nThreads = 8)


net = blockwiseModules(counts, power = 5,
                       TOMType = "unsigned", minModuleSize = 7,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       verbose = 3)


sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
tiff("WGCNA_Dendrogram_23_0831.tiff", res = 300, height = 5.5, width = 7.5, units = "in")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()




moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];




metabolite_metadata <- read.csv("Human_Search/230830_WGCNA_Metadata.csv")

metabolite_metadata$Condition <- NULL
# Make factors numeric
# metabolite_metadata <- data.frame(lapply(metabolite_metadata, as.factor))
# 
# # Convert factor columns to numeric by assigning an integer value to each unique factor level
metabolite_metadata[] <- lapply(metabolite_metadata, function(col) {
  if (is.factor(col)) {
    return(as.numeric(col))
  } else {
    return(col)
  }
})
# 
# # Print the modified data frame
# print(metabolite_metadata)


rownames(metabolite_metadata) <- metabolite_metadata$Name
metabolite_metadata$Name <- NULL


MEs0 = moduleEigengenes(counts, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, metabolite_metadata, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);





sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), " (",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(metabolite_metadata),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

library(ggplot2)
library(dplyr)
library(tidyverse)
meta <- read.csv("Human_Search/230830_WGCNA_Metadata.csv")
meta$Condition <- factor(meta$Condition,levels = c("Cont", "KD", "IR", "IR_KD"))
# Tile Plot of Module-trait Relationships
module_order = names(MEs0) %>% gsub("ME","", .)

MEs0$treatment = meta$Condition
# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  xlab("Subtype") +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
ggsave("23_0831_wgcna_module-trait_relationships.tiff", device = "tiff", 
       dpi = 300, width = 3.75, height = 6)

#Boxplots for each module
for(i in colnames(MEs) ){
  dt <- data.frame(module = MEs[,i], Condition = meta$Condition)
  
  ggplot(dt, aes(x=Condition, y=module, group=Condition, color=Condition)) + 
    geom_boxplot() + 
    scale_x_discrete(labels=c("Cont", "KD", "IR", "IR_KD"))+
    #geom_smooth(formula = time ~ module)
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(-0.3, 1), breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1))
  ggsave(paste0("23_1214_wgcna_module_",i,".tiff"), device = "tiff", 
         dpi = 300, width = 4, height = 4)
}

## Modules of Interest
module_df <- data.frame(
  gene_id = names(net$colors),
  colors = labels2colors(net$colors)
)
## The modules white, lightgreen, green, darkgreen are upregulated
## only in MPBC (lightgreen did not look good in the plot below)
modules_of_interest = c("yellow", "skyblue", "lightcyan", "green", "brown", "black", "red", "pink",
                        "tan", "salmon", "turquoise")

## Pull out list of genes in the modules of interest for MPBC
submod = module_df %>%
  subset(colors %in% modules_of_interest)
row.names(module_df) = module_df$gene_id

subexpr = metabolite_counts[metabolite_counts$PG.UniProtIds %in% submod$gene_id,]
row.names(subexpr) <- subexpr$PG.UniProtIds
subexpr<- subexpr[,-1]
colnames(subexpr) <- meta$Condition

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df$name <- factor(submod_df$name, levels = 
                           c("Cont", "Cont.1", "Cont.2", "Cont.3", "Cont.4", "Cont.5",
                             "Cont.6", "Cont.7", "Cont.8", "Cont.9", "Cont.10", "Cont.11",
                             "Cont.12", "Cont.13", "MPBC", "MPBC.1",  "MPBC.2",  "MPBC.3", 
                             "MPBC.4" , "MPBC.5"  ,"MPBC.6"  ,"MPBC.7"  ,"MPBC.8"  ,"MPBC.9", 
                             "MPBC.10", "MPBC.11", "MPBC.12", "MPBC.13" ,"TNBC"    ,"TNBC.1", 
                             "TNBC.2" , "TNBC.3"  ,"TNBC.4"  ,"TNBC.5"  ,"TNBC.6"  ,"TNBC.7", 
                             "TNBC.8" , "TNBC.9"  ,"TNBC.10" ,"TNBC.11" ,"TNBC.12" ,"TNBC.13",
                             "HER2"   , "HER2.1"  ,"HER2.2"  ,"HER2.3"  ,"HER2.4"  ,"HER2.5", 
                             "HER2.6" , "HER2.7"  ,"HER2.8"  ,"HER2.9"  ,"HER2.10" ,"HER2.11",
                             "HER2.12", "HER2.13","LumB"    ,"LumB.1"  ,"LumB.2"  ,"LumB.3", 
                             "LumB.4" , "LumB.5"  ,"LumB.6"  ,"LumB.7"  ,"LumB.8"  ,"LumB.9", 
                             "LumB.10", "LumB.11", "LumB.12", "LumB.13" ,
                             "LumA"    ,"LumA.1", 
                             "LumA.2" , "LumA.3"  ,"LumA.4"  ,"LumA.5"  ,"LumA.6"  ,"LumA.7", 
                             "LumA.8" , "LumA.9"  ,"LumA.10" ,"LumA.11" ,"LumA.12" ,"LumA.13"
                              ))


submod_df %>% ggplot(., aes(x=name, y=value/1000, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position =  "none"
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "Subtype",
       y = "Abundance (x 10^3)")

ggsave("23_0520_abundance_line_MPBC.tiff", device = "tiff", 
       dpi = 300, width = 3.75, height = 6)
# 
# for(i in unique(submod_df$module)){
#   dt <- submod_df[submod_df$module == i,]
#   
#   dt %>% ggplot(., aes(x=name, y=value/1000, group=gene_id)) +
#     geom_line(aes(color = module),
#               alpha = 0.4) +
#     theme_bw() +
#     theme(
#       axis.text.x = element_text(angle = 90),
#       legend.position =  "none"
#     ) +
#     labs(x = "Subtype",
#          y = "Abundance (x 10^3)")
#   
#   ggsave(paste0("23_0520_abundance_line_",i,".tiff"), device = "tiff", 
#          dpi = 300, width = 3.75, height = 3)
# }
# 




# Export to Cytoscape

genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = metabolite_counts[metabolite_counts$PG.UniProtIds %in% genes_of_interest$gene_id,]
row.names(expr_of_interest) <- expr_of_interest$PG.UniProtIds
expr_of_interest<- expr_of_interest[,-1]
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = 5)

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")



# Heatmaps for each module

# Get proteins for each of the modules of interest
myModule <- cbind(metabolite_counts$PG.UniProtIds, moduleColors)
colnames(myModule) <- c("ID", "moduleColors")
proAssignment <- cbind(Genes = protein$PG.Genes, ID = protein$PG.UniProtIds)
myModule <- data.frame(myModule)
proAssignment <- data.frame(proAssignment)
any(duplicated(myModule$ID))
any(duplicated(proAssignment$ID))
myModule <- merge(proAssignment, myModule, by = "ID")
# Modules of Interest

myModule[moduleColors == "white",1]
myModule[moduleColors == "royalblue",1]
myModule[moduleColors == "green",1]


library(scales)
library(hrbrthemes)
library(tidyr)

# Load the data
all_comparisons <- read.csv("~/Desktop/CRUK_Storming_Cancer/Breast_FFPE/JB8-14_2rep/data/22_1018_JB8-JB14_2rep_panHuman_v01_candidates.csv", stringsAsFactors = F)

# Define the UniProt IDs for each module color
#my_ids <- as.data.frame(myModule) %>%
#  group_by(moduleColors) %>%
#  summarize(UniProtIds = paste(V1, collapse = ",")) %>%
#  mutate(UniProtIds = strsplit(UniProtIds, ","))
#
# Create a plot for each module color
plots_list <- lapply(my_ids$UniProtIds, function(ids) {
  # Filter the data based on UniProt IDs
  cf <- all_comparisons %>%
    filter(UniProtIds %in% ids & Absolute.AVG.Log2.Ratio >= 0.58 & Qvalue <= 0.001) %>%
    mutate(Condition.Numerator = factor(Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR")))
  
  # Create the plot
  p <- ggplot(cf, aes(x = Condition.Numerator, y = Genes, fill = AVG.Log2.Ratio)) +
    coord_equal(ratio = 1/2) +
    geom_tile(color = "white") +
    scale_y_discrete(limits = rev(levels(cf$Genes))) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                         "Log2(FC)", limits = c(-6, 6), oob = squish) +
    #geom_point(aes(x = Condition.Numerator, y = Genes, size = -log10(Qvalue))) +
    #scale_size_binned(range = c(0.5, 3.5)) +
    xlab("Subtype") +
    scale_x_discrete(labels = c("MBC", "TNBC", "HER2", "LumB", "LumA")) +
    ylab("Proteins") +
    theme_ipsum() +
    theme(axis.text = element_text(size = 16, color = "black"), 
          axis.title = element_text(size = 18, color = "black"),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(angle = 90))
  
  # Return the plot
  return(p)
})

# Save the plots
for (i in seq_along(my_ids$UniProtIds)) {
  tiff(paste0("23_0520_Heatmap_", myModule[i, 2], ".tiff"), 
       res = 300, height = 6, width = 3.75, units = "in")
  print(plots_list[[i]])
  dev.off()
}

# Barplots for MPBC proteins of interest

library(matrixStats)
# Barplots for Specific Proteins
progression <- c("Cont", "MPBC", "TNBC", "HER2", "LumB", "LumA")

df <- as.data.frame(my_pro)

mf <- data.frame(Gene = rownames(my_pro),
                 Cont = apply(my_pro[, c(1:14), drop = FALSE], 1, mean),
                 MPBC = apply(my_pro[, c(15:28), drop = FALSE], 1, mean),
                 TNBC = apply(my_pro[, c(71:84), drop = FALSE], 1, mean),
                 HER2 = apply(my_pro[, c(43:56), drop = FALSE], 1, mean),
                 LumB = apply(my_pro[, c(57:70), drop = FALSE], 1, mean),
                 LumA = apply(my_pro[, c(29:42), drop = FALSE], 1, sd)
)

sf <- data.frame(Gene = rownames(my_pro),
                 Cont = apply(my_pro[, c(1:14), drop = FALSE], 1, sd),
                 MPBC = apply(my_pro[, c(15:28), drop = FALSE], 1, sd),
                 TNBC = apply(my_pro[, c(71:84), drop = FALSE], 1, sd),
                 HER2 = apply(my_pro[, c(43:56), drop = FALSE], 1, sd),
                 LumB = apply(my_pro[, c(57:70), drop = FALSE], 1, sd),
                 LumA = apply(my_pro[, c(29:42), drop = FALSE], 1, sd)
)

mf  <-mf %>%  
  pivot_longer(!Gene, names_to = "Group", values_to = "Mean")

sf <- sf %>% pivot_longer(!Gene, names_to = "Group", values_to = "SD")

plotFrame <- full_join(mf, sf)
plotFrame$Group <- factor(plotFrame$Group, levels = progression)



mybar <- function(myGene){
  
  p <- plotFrame %>% filter(Gene == myGene) %>%
    ggplot(aes(x = Group, y = Mean, fill = Group))+
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymax = Mean + SD, ymin = Mean - SD), width = 0.5) +
    ggtitle(myGene) +
    xlab("Subtype") +
    scale_x_discrete(labels = c("Cont", "MPBC", "TNBC", "HER2", "LumB", "LumA"))+
    ylab("Average Abundance") +
    theme_bw()+
    theme(legend.position = "none",axis.text = element_text(size = 16, color = "black"), 
          axis.title = element_text(size = 18, color = "black"),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(angle = 90),
          title = element_text(size = 20, color = "black"))
  
  return(p)
}

# Define the UniProt IDs for each module color
my_ids <- as.data.frame(myModule) %>%
  group_by(moduleColors) %>% filter(moduleColors == "cyan")

plots_list <- lapply(my_ids$Genes, mybar)

unique_values <- unique(my_ids$moduleColors)
for (value in unique_values) {
  folder_name <- paste0("barplots/",value)  # Create folder name with a suffix
  dir.create(folder_name)
}

for (i in seq_along(my_ids$Genes)) {
  tiff(paste0("barplots/",my_ids[i,3], "/","23_0520_barplot_", my_ids[i, 2], ".tiff"), 
       res = 300, height = 3.75, width = 3.75, units = "in")
  print(plots_list[[i]])
  dev.off()
}

#barplots for all points
plotFrame <- data.frame(Gene = rownames(my_pro),
                        my_pro[,c(1:14, 15:28, 71:84, 43:56, 57:70, 29:42)]) %>%
  pivot_longer(-Gene)

plotFrame$name <- factor(plotFrame$name, levels = 
                           c("Cont", "Cont.1", "Cont.2", "Cont.3", "Cont.4", "Cont.5",
                             "Cont.6", "Cont.7", "Cont.8", "Cont.9", "Cont.10", "Cont.11",
                             "Cont.12", "Cont.13", "MPBC", "MPBC.1",  "MPBC.2",  "MPBC.3", 
                             "MPBC.4" , "MPBC.5"  ,"MPBC.6"  ,"MPBC.7"  ,"MPBC.8"  ,"MPBC.9", 
                             "MPBC.10", "MPBC.11", "MPBC.12", "MPBC.13" ,"TNBC"    ,"TNBC.1", 
                             "TNBC.2" , "TNBC.3"  ,"TNBC.4"  ,"TNBC.5"  ,"TNBC.6"  ,"TNBC.7", 
                             "TNBC.8" , "TNBC.9"  ,"TNBC.10" ,"TNBC.11" ,"TNBC.12" ,"TNBC.13",
                             "HER2"   , "HER2.1"  ,"HER2.2"  ,"HER2.3"  ,"HER2.4"  ,"HER2.5", 
                             "HER2.6" , "HER2.7"  ,"HER2.8"  ,"HER2.9"  ,"HER2.10" ,"HER2.11",
                             "HER2.12", "HER2.13","LumB"    ,"LumB.1"  ,"LumB.2"  ,"LumB.3", 
                             "LumB.4" , "LumB.5"  ,"LumB.6"  ,"LumB.7"  ,"LumB.8"  ,"LumB.9", 
                             "LumB.10", "LumB.11", "LumB.12", "LumB.13" ,
                             "LumA"    ,"LumA.1", 
                             "LumA.2" , "LumA.3"  ,"LumA.4"  ,"LumA.5"  ,"LumA.6"  ,"LumA.7", 
                             "LumA.8" , "LumA.9"  ,"LumA.10" ,"LumA.11" ,"LumA.12" ,"LumA.13"
                           ))



mybar <- function(myGene){
  
  p <- plotFrame %>% filter(Gene == myGene) %>%
    ggplot(aes(x = name, y = value))+
    geom_bar(stat = "identity") +
    #geom_errorbar(aes(ymax = Mean + SD, ymin = Mean - SD), width = 0.5) +
    ggtitle(myGene) +
    xlab("Subtype") +
  #  scale_x_discrete(labels = c("Cont", "MPBC", "TNBC", "HER2", "LumB", "LumA"))+
    ylab("Abundance") +
    theme_bw()+
    theme(legend.position = "none",axis.text = element_text(size = 16, color = "black"), 
          axis.title = element_text(size = 18, color = "black"),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(angle = 90),
          title = element_text(size = 20, color = "black"))
  
  return(p)
}

plots_list <- lapply(my_ids$Genes, mybar)

for (i in seq_along(my_ids$Genes)) {
  tiff(paste0("barplots/",my_ids[i,3], "/","23_0520_barplot_", my_ids[i, 2], ".tiff"), 
       res = 300, height = 3.75, width = 3.75, units = "in")
  print(plots_list[[i]])
  dev.off()
}

