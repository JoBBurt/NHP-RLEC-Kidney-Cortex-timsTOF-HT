library(ggplot2)
library(dplyr)
library(scales)
library(hrbrthemes)
library(tidyr)
dq<- all_comparisons

#6.0 SASP
SASP <- c("TNC", "VCAN", "POSTN", "FN1", "FKBP10", "MSN", "HSP90B1", "HSP90AB1", "HSP90AA1", "SERPINH1", "GDI2",
          "VIM", "NID2", "HSPG2", "IGFBP2", "IGFBP4", "KRT2", "KRT10")
cf <- dq %>% filter(Genes %in% SASP)
cf$Condition.Numerator <- factor(cf$Condition.Numerator, levels = c("Meta", "TN", "HER2", "TP", "ERPR"))
cf$Genes <- factor(cf$Genes, levels = SASP)

tiff("output/Heatmap_SASP_2022_0510.tiff", res = 300, height = 7, width = 7, units = "in")
ggplot(cf, aes(x=Condition.Numerator, y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(cf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
 # geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
#  scale_size_binned(range = c(0.25,3))+
  xlab("Subtype") +
  scale_x_discrete(labels = c("MBC", "TNBC", "HER2", "LumB", "LumA")) +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 18, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90))
dev.off()



# Overlap between IR_KD/Cont and others vs Cont
overlap <- c("C4A","CAPN13","CKAP4","CMA1","DES","F13A1","FBLN5","FERMT3","POSTN","SPTA1","TNC","VCAN")
cf <- filter(dq, dq$Genes %in% overlap)
cf$Comparison..group1.group2. <- factor(cf$Comparison..group1.group2., levels = c("IR_KD / Cont", "IR / Cont", "KD / Cont", "IR_KD / IR", "IR_KD / KD", "KD / IR"))

tiff("output/Heatmap_Overlap_2023_0925.tiff", res = 300,units = "in")
ggplot(cf, aes(x=Comparison..group1.group2., y=Genes, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/2) + # aspect ratio y/x
  geom_tile(color = "white") +
  scale_y_discrete(limits = rev(levels(cf$Genes))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       "Log2(FC)",
                       limits=c(-6,6), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  #geom_point(aes(x=Condition.Numerator, y=Genes, size = -log10(Qvalue))) +
  #scale_size_binned(range = c(0.1,3))+
  xlab("Subtype") +
  ylab("Proteins") +
  theme_ipsum() +
  theme(axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))
dev.off()







#modules of interest
dq <- read.csv("output/WGCNA/Homologous_Comparisons_for_WGCNA.csv", stringsAsFactors = F)
modules <- read.csv("output/WGCNA/23_1214_proteins_of_interest_for_modules_of_interest_comparisons.csv", stringsAsFactors = F)
dm <- modules %>% group_by(colors) %>% summarize(Genes = unique(gene_id))
df <- left_join(dm, dq, by = "Genes")
df$Comparison..group1.group2. <- factor(df$Comparison..group1.group2., levels = c("IR_KD / Cont", "IR / Cont", "KD / Cont", "IR_KD / IR", "IR_KD / KD", "KD / IR"))
dff <- df %>% filter(!(AVG.Log2.Ratio >= -0.58 & AVG.Log2.Ratio <= 0.58 & Qvalue > 0.01))
# Define the set of values you want to keep
values_to_keep <- c("IR_KD / IR", "KD / Cont", "IR / Cont")  # Replace with your actual values

# Sample condition function: modify according to your specific condition
condition_function <- function(x) {
  any(x > threshold)  # Replace with your actual condition
}
threshold <- 0.58
# Applying the filter
filtered_df <- df %>%
  filter(Comparison..group1.group2. %in% values_to_keep) %>%
  group_by(Genes) %>%
  mutate(condition_met = condition_function(Absolute.AVG.Log2.Ratio)) %>%
  filter(condition_met) %>%
  ungroup()  # Removing the grouping


# Loop through unique colors
for(i in unique(dm$colors)){
  print(paste("Processing color:", i))
  
  # Filter the DataFrame
  cf <- filtered_df %>% filter(colors == i)
  
  # Debug: Print head of cf
  print(head(cf))
  
  # Construct the TIFF file name
  tiff_filename <- paste("output/WGCNA/Heatmap_", i, "_filtered_2023_1214.tiff", sep = "")
  print(paste("Saving to:", tiff_filename))
  
  # Generate and save the plot
  #tiff(tiff_filename, res = 300, height = 4, width = 4, units = "in")
  print("Creating plot...")
 p <- cf %>%
    ggplot(., aes(x=Comparison..group1.group2., y=Genes, fill=AVG.Log2.Ratio)) +
    coord_fixed(ratio = 1) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                         "Log2(FC)",
                         limits=c(-3.5,3.5), oob=squish) +
    xlab("Comparison") +
    scale_x_discrete(labels = c("IR_KD / Cont", "IR / Cont", "KD / Cont", "IR_KD / IR", "IR_KD / KD", "KD / IR")) +
    ylab("Proteins") +
    theme_ipsum() +
    theme(axis.text = element_text(size = 14, color = "black"), 
          axis.title = element_text(size = 18, color = "black"),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(angle = 90))
  #dev.off()
  
  # Save the plot using ggsave
  ggsave(filename = tiff_filename, plot = p, dpi = 300, height = 22, width = 8, units = "in")
  
  print(paste("Plot saved as:", tiff_filename))
  print("Plot saved.")
}



