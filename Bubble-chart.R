#Making Proportion of superclasses per apple figure

library(ggplot2)
library(reshape2) 

#Generate table consisting of columns for x-axis, y-axis and proportions:

##TREATMENT OF View.LABORATORY FUNGI
feature.to.fungus <- filter(data.filt.blanks, ATTRIBUTE_Type == "Laboratory")
fungus.features <- filter(feature.to.fungus) %>%
  select(-ATTRIBUTE_Genus)
fungus.features <- fungus.features[,c(7:ncol(fungus.features))]
fungus.features <- t(fungus.features)
fungus.features <- as.data.frame(fungus.features)
fungus.features$row.ID2 <- rownames(fungus.features)
fungus.features <- merge(data3, fungus.features, by = "row.ID2")
fungus.features <- select(fungus.features, -row.ID2)
fungus.features$rowsums <- rowSums(fungus.features[2:ncol(fungus.features)])

features.detected.fungus <- select(fungus.features, feature.ID, rowsums) %>%
  filter(rowsums > 0)
featureIDs <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/feature_id_to_mass_pre-tidied.csv")
features.detected.fungus <- merge(features.detected.fungus, featureIDs, by = "feature.ID")

##TREATMENT OF APPLES
apple <- filter(data.filt.blanks, ATTRIBUTE_Type == "Apple_natural" | ATTRIBUTE_Type == "Apple_excision")
apple_target_list <- apple$filename

apple.superclasses <- function (apples) {
apple <- filter(data.filt.blanks, ATTRIBUTE_Type == "Apple_natural" | ATTRIBUTE_Type == "Apple_excision")
apple <- filter(apple, filename == as.character(apples)) ## exclude if not looking at specific sample

apple.feat.only <- apple[,-c(1:7)]
apple.features <- t(apple.feat.only)
apple.features <- as.data.frame(apple.features)
apple.features$row.ID2 <- rownames(apple.features)

apple.features <- merge(data3, apple.features, by = "row.ID2") %>%
  filter(V1 > 0)
apple.features <- select(apple.features, -row.ID2)

features.detected.apples <- select(apple.features, feature.ID)
features.detected <- merge(features.detected.fungus, features.detected.apples, by = "feature.ID")

in.both2 <- merge(app.feats, features.detected, by = "feature.ID")


df <- in.both2 %>% 
  group_by(NPC_NPC.superclass) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  filter(!is.na(NPC_NPC.superclass)) %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(round(perc, 3))) %>%
  mutate(lab.ypos = cumsum(perc) - 0.5*perc)
df$x <- paste0(df$NPC_NPC.superclass, " (", df$labels, ")")
print(df)
return(df)
}

#Calculates the superclasses of each individual apple
df3 <- data.frame(matrix(ncol=3, nrow = 0)) ##sets up a new df
for (i in apple_target_list) {
  print(i)
  df <- apple.superclasses(i)
  df2 <- data.frame(matrix(ncol=0, nrow = nrow(df))) ##sets up a new df
  if (nrow(df) > 0) {
    df2 <- cbind(df2, df$NPC_NPC.superclass, df$perc, i)
  }
  df3 <- rbind(df3, df2)
}

colnames(df3)[3] <- "filename"

appleIDs <- read.csv("appleIDs.csv")
merged_df <- merge(df3, appleIDs, by = "filename") #converts filename to appleIDs

#Generate bubble chart figure

category_colors <- rainbow(length(unique(merged_df$`df$NPC_NPC.superclass`)))

plot_obj  <- ggplot(merged_df, aes(x = `appleID`, y = `df$NPC_NPC.superclass`, size = `df$perc`, color = `df$NPC_NPC.superclass`)) +
  geom_point() +
  scale_size(range = c(1, 8)) +
  scale_color_manual(values = category_colors) +  # Set the color scale
  labs(x = "Black Apple ID", y = "Superclasses")

plot_obj + guides(color = FALSE) + labs(size = "Proportion of Metabolites") + 
  theme_minimal() + theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_y_discrete(limits = rev(levels(plot_obj$data$`df$NPC_NPC.superclass`)))
#Change width to 1500