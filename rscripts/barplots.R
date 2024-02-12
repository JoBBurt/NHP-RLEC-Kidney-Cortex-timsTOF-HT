library(matrixStats)
# Barplots for Specific Proteins
progression <- c("Hela", "HeLaFlagRECs1", "Chloro","Doxo", "D+C")

df <- as.data.frame(my_pro)

mf <- sapply(split.default(df, names(df)), rowMeans)
mf <- mf[,progression]
mf <- as.data.frame(mf)
sf <- data.frame(Gene = rownames(my_pro),
                 DF = rowSds(my_pro[,c(rep("Hela", 3))]),
                 Meta = rowSds(my_pro[,c(rep("HeLaFlagRECs1", 3))]),
                 TN = rowSds(my_pro[,c(rep("Chloro", 3))]),
                 HER2 = rowSds(my_pro[,c(rep("Doxo", 3))]),
                 ERPR = rowSds(my_pro[,c(rep("D+C", 3))])
)

mf  <-mf %>%  mutate(Gene = rownames(mf)) %>% 
  pivot_longer(!Gene, names_to = "Group", values_to = "Mean")

sf <- sf %>% pivot_longer(!Gene, names_to = "Group", values_to = "SD")

plotFrame <- full_join(mf, sf)
plotFrame$Group <- factor(plotFrame$Group, levels = progression)


mybar <- function(myGene){
    
  plotFrame %>% filter(Gene == myGene, Group != "NA") %>%
    ggplot(aes(x = Group, y = Mean, fill = Group))+
    geom_bar(stat = "identity", fill = mycol) +
    geom_errorbar(aes(ymax = Mean + SD, ymin = Mean - SD), width = 1) +
    ggtitle(myGene) +
    xlab("Treatment") +
    scale_x_discrete(labels = c("WT", "RECs1", "Chloro", "Doxo", "D+C"))+
    ylab("Average Abundance") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          axis.text = element_text(size = 16, color = "black"), 
          axis.title = element_text(size = 18, color = "black"),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(angle = 90),
          title = element_text(size = 20, color = "black"))
  
  ggsave(filename = paste(myGene, "Barplot.tiff"), 
         device = tiff, path = "output", dpi = "print", width = 5, height = 5, units = 'in')  
  
}

mycol <- c("red", "darkgreen", "blue", "orange", "purple")
mybar("SERPINH1")
mybar("PRELP")
mybar("POSTN")
mybar("HSPG2")
mybar("TIMP1")
mybar("GPC1")
mybar("TGFB1")
mybar("SELENBP1")
mybar("TMBIM1")
