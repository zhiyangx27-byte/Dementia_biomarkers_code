library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

omics_data_df <- read.csv("/Users/a1234/Desktop/omics_data(class&sub).csv")

dementia_data <- omics_data_df %>%
  filter(Cognitive.Status == "No dementia")

class_proportion <- dementia_data %>%
  group_by(Class) %>%
  summarise(count = n()) %>%
  mutate(ratio = count / sum(count)) 

num_classes <- length(unique(class_proportion$Class))
color_palette_class <- colorRampPalette(brewer.pal(12, "Set3"))(num_classes)
names(color_palette_class) <- unique(class_proportion$Class)

class_bar <- ggplot(class_proportion, aes(x = Class, y = ratio, fill = Class)) +
  geom_col(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = color_palette_class) +
  labs(title = "Proportion of Each Class in Dementia Cells",
       x = "Class",
       y = "Ratio") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  guides(fill = guide_legend(ncol = 1))

print(class_bar)

threshold <- 0.02
subclass_proportion <- dementia_data %>%
  group_by(Subclass) %>%
  summarise(count = n()) %>%
  mutate(ratio = count / sum(count)) %>%
  mutate(Subclass = as.character(Subclass)) %>%
  mutate(Subclass = ifelse(ratio < threshold, "Others(< 2.00%)", Subclass)) %>%
  group_by(Subclass) %>%
  summarise(ratio = sum(ratio), .groups = "drop") 

color_palette_subclass <- colorRampPalette(brewer.pal(12, "Set3"))(13)
names(color_palette_subclass) <- levels(factor(subclass_proportion$Subclass))

subclass_bar <- ggplot(subclass_proportion, aes(x = Subclass, y = ratio, fill = Subclass)) +
  geom_col(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = color_palette_subclass) +
  labs(#title = "Subclass Distribution in Dementia Cells",
       x = "Subclass",
       y = "Ratio") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  guides(fill = guide_legend(ncol = 1))

print(subclass_bar)

