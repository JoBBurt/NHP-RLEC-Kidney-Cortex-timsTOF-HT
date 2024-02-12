## Set the working directory
setwd("~/Desktop/RCCN_Kidney/Cortex/Batch_1and2/")

## Load search data from Macaca mullata proteome

## This is the treatment factor used in the PCA and heatmap
Rhe_Con <- read.csv("Rhesus_Search/230602_RCCN2_bothAB_Kidney_Cortex_directDIA_v8_ConditionSetup_v02.csv",
                      stringsAsFactors = F)
Rhe_Pro <- read.csv("Rhesus_Search/20230612_RCCN2_bothAB_Kidney_Cortex_directDIA_v8_Protein_Report_2pep.csv",
                    stringsAsFactors = F)
Rhe_Pep <- read.csv("Rhesus_Search/20230719_131725_230602_RCCN2_bothAB_Kidney_Cortex_directDIA_v8_Peptide_Quant_Report_2pep.csv",
                    stringsAsFactors = F)
Rhe_Pair <- read.csv("Rhesus_Search/230612_RCCN2_bothAB_Kidney_Cortex_directDIA_v8_candidates_2pep_V03.csv",
                            stringsAsFactors = F)

## Load search data from Human proteome

Hum_Con <- read.csv("Human_Search/230602_RCCN2_bothAB_Kidney_Cortex_ConditionSetup_v02.csv",
                             stringsAsFactors = F)
Hum_Pro <- read.csv("Human_Search/20230718_171725_230602_RCCN2_bothAB_Kidney_Cortex_HumanFASTA_directDIA_v9_Protien_Quant_2pep_Report.csv",
                    stringsAsFactors = F)
Hum_Pep <- read.csv("Human_Search/20230718_172014_230602_RCCN2_bothAB_Kidney_Cortex_HumanFASTA_directDIA_v9__Peptide_Quant_2pep_Report.csv",
                    stringsAsFactors = F)
Hum_Pair <- read.csv("Human_Search/230612_RCCN2_bothAB_Kidney_Cortex_HumanFASTA_directDIA_v9_candidates_2pep_v03.csv",
                            stringsAsFactors = F)
# Join data together
library(dplyr)
Homolog_Pep <- inner_join(Rhe_Pep, Hum_Pep, by = "EG.PrecursorId")
#Homolog_Pro <- inner_join(Rhe_Pro, Hum_Pro, by = "PG.ProteinDescriptions")
temp <- unique(Homolog_Pep$PG.Genes.y) 
Homolog_Pro2 <- left_join(Rhe_Pro, Homolog_Pep, by = c("PG.Genes" = "PG.Genes.y"))
Homolog_Pro3 <- left_join(Homolog_Pro2, Homolog_Pep, by = c("PG.Genes" = "PG.Genes.x"), relationship = "many-to-many")
Homolog_Pro3 <- Homolog_Pro3[, c(1:8, 390:396)]
temp <- unique(Homolog_Pro3$PG.Genes) 

write.csv(Homolog_Pep, "output/Data_Tables/Homolog_Pep_V01_1pep.csv")
write.csv(Homolog_Pro3, "output/Data_Tables/Homolog_Pro_V02.csv")

head(Hum_Pep$EG.PrecursorId)

#1 pep

Rhe_Pep <- read.csv("Rhesus_Search/230602_RCCN2_bothAB_Kidney_Cortex_directDIA_v8_Report_Birgit_Peptide_Quant_1pep.csv",
                    stringsAsFactors = F)
Rhe_Pro <- read.csv("Rhesus_Search/230602_RCCN2_bothAB_Kidney_Cortex_directDIA_v8_Report_Birgit_Protein_Quant_Pivot_1pep.csv",
                    stringsAsFactors = F)
Homolog_Pep <- inner_join(Rhe_Pep, Hum_Pep, by = "EG.PrecursorId")
#Homolog_Pro <- inner_join(Rhe_Pro, Hum_Pro, by = "PG.ProteinDescriptions")
temp <- unique(Homolog_Pep$PG.Genes.y) 
Homolog_Pro2 <- left_join(Rhe_Pro, Homolog_Pep, by = c("PG.Genes" = "PG.Genes.y"))
Homolog_Pro3 <- left_join(Homolog_Pro2, Homolog_Pep, by = c("PG.Genes" = "PG.Genes.x"), relationship = "many-to-many")
Homolog_Pro3 <- Homolog_Pro3[, c(1:8, 390:396)]
temp <- unique(Homolog_Pro3$PG.Genes) 
write.csv(Homolog_Pro3, "output/Data_Tables/Homolog_Pro_1pep_V02.csv")


pep_1 <- read.csv("output/Data_Tables/Homolog_Pro_1pep_v02.csv", stringsAsFactors = F)
pep_2 <- read.csv("output/Data_Tables/Homolog_raw/Homolog_Pro_v02.csv", stringsAsFactors = F)

temp <- anti_join(pep_1, pep_2, by = "PG.Genes.y")

write.csv(temp, "output/Data_Tables/Holog_Pro_1pep_V03.csv")


temp <- read.csv("output/temp.csv", stringsAsFactors = F)
head(temp)
non<- anti_join(Rhe_Pro, temp, by = "PG.Genes")


#
# I need to filter out the proteins that we did not find in the overlap between human and 
# NHPs 

# Use the human proteins in homology excel sheet as the key for human comparison file 
# and protein file

key <- read.csv("output/Data_Tables/homology_key.csv", stringsAsFactors = F)
all_comparisons <- read.csv("Human_Search/230612_RCCN2_bothAB_Kidney_Cortex_HumanFASTA_directDIA_v9_candidates_2pep_v03.csv",
                            stringsAsFactors = F)


homo_comparisons <- all_comparisons %>%
  filter(Genes %in% key$HS.Genes)


# Make a histogram of peptides/protein that were matched between the two species
Homolog_Pep$PeptideID <- sapply(strsplit(Homolog_Pep$EG.PrecursorId, "\\."), '[', 1) 

NHP_count <- Homolog_Pep %>% group_by(PG.Genes.x) %>% summarize(Count = n())
Human_count <- Homolog_Pep %>% group_by(PG.Genes.y) %>% summarize(Count = n())

ggplot(NHP_count, aes(x = Count)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of NHP Peptide IDs per NHP Gene",
       x = "Number of Peptide IDs",
       y = "Frequency") +
  theme_minimal()


ggplot(Human_count, aes(x = Count)) +
  geom_histogram(binwidth = 1, fill = "#E9FFE9", color = "black") +
  labs(title = "Histogram of NHP Peptide IDs per Human Gene",
       x = "Number of Peptide IDs",
       y = "Frequency") +
  xlim(c(0,100))+
  #facet_wrap(~factor(Comparison..group1.group2.), nrow = 2) +
  theme_classic() +
  theme(axis.title =element_text(size = 20, color = "black"),
        axis.text = element_text(size = 18, color = "black"),
        legend.position = "none",
        title = element_blank())
ggsave(filename = "NHP_Peptides_per_Human_Gene.tiff",
        device = tiff, path = "output", dpi = "print", width = 5, height = 5, units = 'in') 
