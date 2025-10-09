library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

omics_data_df <- read.csv("/Users/a1234/Desktop/omics_data(class&sub).csv")

dementia_data <- omics_data_df %>%
  filter(Cognitive.Status == "Dementia")

class_proportion <- dementia_data %>%
  group_by(Class) %>%
  summarise(count = n()) %>%
  mutate(ratio = count / sum(count)) %>%
  mutate(label = paste0(Class, "\n(", round(ratio * 100, 2), "%)\n(", count, ")"))

num_classes <- length(unique(class_proportion$Class))
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(num_classes)
names(color_palette) <- unique(class_proportion$Class)

class_pie <- ggplot(class_proportion, aes(x = "", y = ratio, fill = Class)) +
  geom_col(width = 1, color = "white", linewidth = 0.5) + 
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_palette) +
  geom_label_repel(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.5,  
    min.segment.length = 0.2,  
    box.padding = 0.5,
    max.overlaps = 30,  
    segment.color = "gray40", 
    show.legend = FALSE
  ) +
  theme_void() +
  theme(
    legend.position = "right",  
    legend.text = element_text(size = 8),  
    plot.title = element_text(hjust = 0.5, size = 14)  
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(title = "Proportion of Each Class in Dementia Cells")

print(class_pie)

threshold <- 0.02
subclass_proportion <- dementia_data %>%
  group_by(Subclass) %>%
  summarise(count = n()) %>%
  mutate(ratio = count / sum(count)) %>%
  mutate(Subclass = as.character(Subclass)) %>%
  mutate(Subclass = ifelse(ratio < threshold, "Others(< 2.00%)", Subclass)) %>%
  group_by(Subclass) %>%
  summarise(ratio = sum(ratio), .groups = "drop") %>%
  mutate(label = paste0(Subclass, "\n", scales::percent(ratio)))

color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(13)
names(color_palette) <- levels(factor(subclass_proportion$Subclass))

subclass_pie <- ggplot(subclass_proportion, aes(x = "", y = ratio, fill = Subclass)) +
  geom_col(width = 1, color = "white", linewidth = 0.5) + 
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_palette) +
  geom_label_repel(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    min.segment.length = 0.2,
    box.padding = 0.5,
    max.overlaps = 30,
    segment.color = "gray40",
    show.legend = FALSE
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(title = "Subclass Distribution in Non-dementia Cells")

print(subclass_pie)


