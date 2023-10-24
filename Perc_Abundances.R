# Normalize each row to its maximum value

data.filt.norm <- data.filt.blanks # Make a copy of the original data frame
cols_to_normalize <- 8:ncol(data.filt.norm) # Define the columns to normalize

data.filt.norm[, cols_to_normalize] <- t(apply(data.filt.norm[, cols_to_normalize], 1, function(row) {
  max_val <- max(row, na.rm = TRUE)
  if (max_val != 0) {
    row / max_val
  } else {
    row
  }
}))

#Now to take averages of each row where a metabolite is present, i.e. ignoring zero values

apple <- filter(data.filt.norm, ATTRIBUTE_Type == "Apple_natural" | ATTRIBUTE_Type == "Apple_excision") #Can edit if only interested in excisions of whole black apples
apple.feat.only <- apple[,-c(1:7)]
apple.feat.only[apple.feat.only < 0.0000001] <- NA # Replace zero values with NA in apple.feat.only
apple.averages <- as.data.frame(colMeans(apple.feat.only, na.rm = TRUE)) # Calculate column means while ignoring NA values
row.ID2 <- rownames(apple.averages) ##not comparable with the above data3 object
apple.averages <- cbind(row.ID2,apple.averages)

#Repeating for the lab

lab <- filter(data.filt.norm, ATTRIBUTE_Type == "Laboratory") 
lab.feat.only <- lab[,-c(1:7)]
lab.feat.only[lab.feat.only < 0.0000001] <- NA 
lab.averages <- as.data.frame(colMeans(lab.feat.only, na.rm = TRUE))
row.ID2 <- rownames(lab.averages) ##not comparable with the above data3 object
lab.averages <- cbind(row.ID2,lab.averages)

in.both <- merge(data3, apple.averages, by = "row.ID2")
in.both <- merge(in.both, lab.averages, by = "row.ID2")
# Rename column 3 to "apple.avg" and column 4 to "lab.avg"
colnames(in.both)[3] <- "apple.avg"
colnames(in.both)[4] <- "lab.avg"

#Make ratio of apple to lab
in.both <- mutate(in.both, ratio = apple.avg / lab.avg)
merged.df <- merge(featureIDs, in.both, by = "feature.ID") ## now feature.ID is compatible.
merged.df.with.metadata2 <- merge(app.feats, merged.df, by = "feature.ID")

#Now last column will show the ratio of apple/lab. 