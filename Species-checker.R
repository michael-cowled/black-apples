##FOR UNTARGETED ANALYSIS WITH OSMAC COLLECTION OF METABOLITES
#To find a target list for each species based on frequency of feature in samples (to avoid potential carry-over)

species_data <- data.filt.blanks[,-c(1,2,4:7)] %>%
  filter(ATTRIBUTE_Species == "polonicum") %>% ##change to species for function
  select(-ATTRIBUTE_Species) %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  colSums() %>%
  as.data.frame() %>%
  filter(. >= 3) ##filters for features that are present in at least 3 samples
species_target_list <- row.names(species_data)  

##Changed above code to work with row_ID instead of feature ID

apple_data <- apple[,-c(2:7)]

apple.finder <- function (features) {
  apple.IDs <- select(apple_data, (as.character(features)))
  return(apple.IDs)
}

df <- data.frame(matrix(ncol=0, nrow = 43)) ##sets up a new df
for (i in species_target_list) {
  df <- cbind(df, apple.finder(i))
}

#Calculates the number of apples with selected features:
df <- mutate_if(df, is.numeric, ~1 * (. > 0))
df <- mutate(df, Total = rowSums(df))
paste0(sum(df$Total >= 6), " of total ", nrow(df), " apples") ##Need at least 5 features to be considered present
pol <- df$Total >= 6


##Heatmap of which fungus is present in which apple
#export as h:166, w.800
black.apples <- read.csv("Apples with Fungus present.csv")
ggplot(black.apples, aes(x, y, fill= z)) + 
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  scale_fill_gradient(low="white", high="black") +
  labs(x = "Black Apple ID", y = "Species") +
  theme(legend.position = "none", axis.text.y = element_text(face = "italic"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##Barchart to show proportion/frequency
black.apples2 <- read.csv("Apples with Fungus present2.csv")

p<-ggplot(data=black.apples2, aes(x=x, y=y)) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Species") +
  ylab("Frequency") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"))
p