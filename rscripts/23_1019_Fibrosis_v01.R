library(ggplot2)
library(dplyr)
library(scales)
library(hrbrthemes)
library(tidyr)

proteins <- read.csv("data/23_1019_homologous_protien_quant_v01.csv")
comparisons <- read.csv("data/23_1019_homolog_comparisons_v01.csv")
metadata <- read.csv("data/230602_RCCN2_bothAB_Kidney_Cortex_ConditionSetup_v03.csv")

mypattern <- "RCCN\\d-\\d+"
## my_pro is used for the PCA and Heatmap
my_pro <- prepmat(proteins, mypattern)
colnames(my_pro) <- metadata$ID
rownames(my_pro) <- proteins$PG.Genes

df <- subset(my_pro, rownames(my_pro) %in% overlap)
# Transpose the data frame
transposed_df <- as.data.frame(t(as.matrix(df)))

# Optionally, you might want to update row and column names
rownames(transposed_df) <- colnames(df)
colnames(transposed_df) <- rownames(df)

transposed_df %>%
ggplot(aes(rownames(transposed_df), F13A1, color = metadata$Condition, shape = factor(metadata$Fib))) +
  geom_point() +
  xlab(label = metadata$Condition) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

