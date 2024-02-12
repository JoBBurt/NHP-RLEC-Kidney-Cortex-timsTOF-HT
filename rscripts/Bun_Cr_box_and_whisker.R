# Load necessary libraries
library(ggplot2)
library(readr)

# Load the data
data <- read_csv('Buck_Kidney_Clinical_Measures.csv')

# Create the box and whisker plot
ggplot(data, aes(x = `Kidney Group`, y = BUN)) +
  geom_boxplot(aes(fill = `Kidney Group`)) +
  geom_jitter(color = "black", size = 0.4, width = 0.2) +
  theme_bw() +
  theme(legend.position = 'none')+
  labs(x = "Kidney Group",
       y = "BUN")

# Save the plot
ggsave("BUN_Levels_Boxplot.png", width = 3, height = 3.5)

ggplot(data, aes(x = `Kidney Group`, y = Cr)) +
  geom_boxplot(aes(fill = `Kidney Group`)) +
  geom_jitter(color = "black", size = 0.4, width = 0.2) +
  theme_bw() +
  theme(legend.position = 'none')+
  labs(x = "Kidney Group",
       y = "Creatinine") + 
  ylim(c(0,3))

# Save the plot
ggsave("Cr_Levels_Boxplot.png", width = 3, height = 3.5)
