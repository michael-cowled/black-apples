##FOR TARGETED ANALYSIS FOR KNOWN FEATURES
#A script to identify a fungal species present in apples:

library(dplyr)

apple_data <- apple[,-c(2:7)]

##Convert targeted feature.ID to row.ID2

feature.finder <- function (ID) {
  filtered.data3 <- filter(data3, feature.ID == ID)
  row_no <- filtered.data3$row.ID2
  return(row_no)
}

##Change features to ones of interest!
list.of.features <- c(9126, 9053, 6614, 6975)

apple.finder <- function (features) {
  row_no <- feature.finder(features) ##change to selected features
  apple.IDs <- select(apple_data, (as.character(row_no)))
  return(apple.IDs)
}

df <- data.frame(matrix(ncol=0, nrow = 43)) ##sets up a new df
for (i in list.of.features) {
  df <- cbind(df, apple.finder(i))
}

#Calculates the number of apples with selected features:
colSums(mutate(df, Total = rowSums(df) >0))