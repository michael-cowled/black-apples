library(dplyr)
library(tidyr)
library(stringr)
library(jsonlite)

#Importing NP databases
NPatlas <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/NPdatabases/NPatlas.csv")
NPatlas <- mutate(NPatlas, across('compound_smiles', str_replace_all, '@', ''))
colnames(NPatlas)[c(3, 11)] <- c("compound", "smiles_nostereo")

##EXTRACTION OF DATA FROM MIBIG DATABASE

#Reading and merging mibig compounds
filenames <- list.files(path = "C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/NPdatabases/json/", 
                        pattern = "", full.names = TRUE)
df <- data.frame() #arbitrary df setup
i <- 1 # Initialize the index variable i to start at 1

# Loop through the filenames until i reaches the length of the 'filenames' vector
while (i <= length(filenames)) {
  # Read the JSON file at the current index i and convert it to a data frame with flatten and simplifyVector options
  tempdf <- fromJSON(filenames[i], flatten = TRUE, simplifyVector = TRUE)
  # Extract the 'cluster' component from the data frame and assign it to 'tempdf'
  tempdf <- tempdf$cluster
  # Extract the 'compounds' component from the data frame and assign it to 'tempdf'
  tempdf <- tempdf$compounds
  df <- bind_rows(df, tempdf)   # Bind the rows of 'tempdf' to the data frame 'df'
  i <- i + 1 # Increment the index 'i' by 1 for the next iteration
  # Print the current value of 'i'
}

mibig <- df #renames df of collated compounds from mibig
#Removal of stereochemistry from SMILES
mibig <- mutate(mibig, across('chem_struct', str_replace_all, '@', ''))
colnames(mibig)[3] <- "smiles_nostereo"

#Importing data
#note data file is already pre-tidied with the relatively-empty columns deleted before importing
#m/z and retention time copied to a new table titled "feature_id_to_mass" with the corresponding "ID"
#Make sure each annotation file is in .csv format and has a column name "feature.ID"

metadata <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/metadata.csv")
featureIDs <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/feature_id_to_mass_pre-tidied.csv")
gnpsIDs <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/Annotations/GNPS_annotations.csv")
dereplicatorIDs <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/Annotations/Dereplicator_annotations.csv")
pcdlIDs <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/Annotations/PCDL_annotations.csv")
CSIIDs <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/Annotations/CSI_finger_IDs.csv")

#Tidy CSIIDs
CSIIDs <- filter(CSIIDs, ConfidenceScore >= 0.95 & ConfidenceScore != "N/A")
CSIIDs <- separate(CSIIDs, col = id, into = c("siriusID", "sirius", "feature.ID"), sep = "\\_")
CSIIDs <- CSIIDs[, c(5:9, 11:14, 23)]
colnames(CSIIDs) <- paste("CSI", colnames(CSIIDs), sep = "_")
colnames(CSIIDs)[10] <- "feature.ID"

#Tidy gnpsIDs
gnpsIDs <- gnpsIDs[,c(2,8,9,14, 21, 41:45)]
colnames(gnpsIDs) <- paste("GNPS", colnames(gnpsIDs), sep = "_")
colnames(gnpsIDs)[5] <- "feature.ID"

#Tidy dereplicatorIDs
dereplicatorIDs <- filter(dereplicatorIDs, FDR <=5) #Filters by FDR
dereplicatorIDs <- dereplicatorIDs[,c(4, 6, 14)]
colnames(dereplicatorIDs) <- paste("DEREP", colnames(dereplicatorIDs), sep = "_")
colnames(dereplicatorIDs)[1] <- "feature.ID"

#Tidy pcdlIDs
pcdlIDs <- pcdlIDs[,c(1, 5, 6)]
colnames(pcdlIDs) <- paste("PCDL", colnames(pcdlIDs), sep = "_")
colnames(pcdlIDs)[1] <- "feature.ID"

#Merge all tables
data.merged <- merge(featureIDs, gnpsIDs,
                     by = 'feature.ID', all = TRUE)
data.merged <- merge(data.merged, dereplicatorIDs,
                     by = 'feature.ID', all = TRUE)
data.merged <- merge(data.merged, pcdlIDs,
                     by = 'feature.ID', all = TRUE)
data.merged <- merge(data.merged, CSIIDs,
                     by = 'feature.ID', all = TRUE)
data.merged <- mutate(data.merged, across('GNPS_Smiles', str_replace_all, '@', ''))
data.merged <- mutate(data.merged, across('DEREP_SMILES', str_replace_all, '@', ''))
colnames(data.merged)[c(8, 15)] <- c("GNPS_smiles_nostereo", "DEREP_smiles_nostereo")

##Filter for all features with an annotation
data.merged.filtered <- filter(data.merged, !is.na(GNPS_Compound_Name) | !is.na(GNPS_smiles_nostereo)
                               | !is.na(DEREP_Name) | !is.na(DEREP_smiles_nostereo)
                               | !is.na(PCDL_Compound_Name) | !is.na(CSI_name)
                                        | !is.na(CSI_smiles))

#Fix compound names
n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("Spectral match", data.merged.filtered$GNPS_Compound_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- gsub("Spectral Match to ","", data.merged.filtered$GNPS_Compound_Name[gnps.peak])
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- gsub(" from NIST14","", data.merged.filtered$GNPS_Compound_Name[gnps.peak])
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("Massbank", data.merged.filtered$GNPS_Compound_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- gsub("^.{0,18}","", data.merged.filtered$GNPS_Compound_Name[gnps.peak])
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("NCGC", data.merged.filtered$GNPS_Compound_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- gsub("^.{0,16}","", data.merged.filtered$GNPS_Compound_Name[gnps.peak])
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("40.0 eV", data.merged.filtered$GNPS_Compound_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- substring(data.merged.filtered$GNPS_Compound_Name[gnps.peak], 1, nchar(data.merged.filtered$GNPS_Compound_Name[gnps.peak])-10)
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl(".00 eV", data.merged.filtered$GNPS_Compound_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- substring(data.merged.filtered$GNPS_Compound_Name[gnps.peak], 1, nchar(data.merged.filtered$GNPS_Compound_Name[gnps.peak])-11)
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl('"', data.merged.filtered$GNPS_Compound_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- gsub('"','', data.merged.filtered$GNPS_Compound_Name[gnps.peak])
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("CollisionEnergy", data.merged.filtered$GNPS_Compound_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- substring(data.merged.filtered$GNPS_Compound_Name[gnps.peak], 1, nchar(data.merged.filtered$GNPS_Compound_Name[gnps.peak])-23)
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl(".beta.-", data.merged.filtered$GNPS_Compound_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- gsub(".beta.-","", data.merged.filtered$GNPS_Compound_Name[gnps.peak])
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("|", data.merged.filtered$GNPS_Compound_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- gsub("\\|.*","", data.merged.filtered$GNPS_Compound_Name[gnps.peak])
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("ReSpect", data.merged.filtered$GNPS_Compound_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- gsub("^.{0,17}","", data.merged.filtered$GNPS_Compound_Name[gnps.peak])
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("_", data.merged.filtered$DEREP_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$DEREP_Name[gnps.peak] <- gsub("_"," ", data.merged.filtered$DEREP_Name[gnps.peak])
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("]", data.merged.filtered$GNPS_Compound_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- gsub("\\[", "(", gsub("\\]", ")", data.merged.filtered$GNPS_Compound_Name[gnps.peak]))
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("}", data.merged.filtered$GNPS_Compound_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_Compound_Name[gnps.peak] <- gsub("\\{", "(", gsub("\\}", ")", data.merged.filtered$GNPS_Compound_Name[gnps.peak]))
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("]", data.merged.filtered$CSI_name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$CSI_name[gnps.peak] <- gsub("\\[", "(", gsub("\\]", ")", data.merged.filtered$CSI_name[gnps.peak]))
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("}", data.merged.filtered$CSI_name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$CSI_name[gnps.peak] <- gsub("\\{", "(", gsub("\\}", ")", data.merged.filtered$CSI_name[gnps.peak]))
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("]", data.merged.filtered$DEREP_Name[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$DEREP_Name[gnps.peak] <- gsub("\\[", "(", gsub("\\]", ")", data.merged.filtered$DEREP_Name[gnps.peak]))
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("]", data.merged.filtered$DEREP_smiles_nostereo[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$DEREP_smiles_nostereo[gnps.peak] <- gsub("\\[", "(", gsub("\\]", ")", data.merged.filtered$DEREP_smiles_nostereo[gnps.peak]))
  }
  gnps.peak <- gnps.peak + 1
}

n <- nrow(data.merged.filtered)
gnps.peak <- 1
while (gnps.peak <= n) {
  if ((any(grepl("+1", data.merged.filtered$GNPS_smiles_nostereo[gnps.peak], ignore.case = TRUE), na.rm = TRUE))) {
    data.merged.filtered$GNPS_smiles_nostereo[gnps.peak] <- gsub("+1", "", gsub("-1", "", data.merged.filtered$GNPS_smiles_nostereo[gnps.peak]))
  }
  gnps.peak <- gnps.peak + 1
}

##Filter for all features with an annotation corresponding to a database
#check gnps

df <- data.frame() #arbitrary df setup

i <- 1
while (i <= nrow(data.merged.filtered)) {

if (any(grepl(data.merged.filtered$GNPS_Compound_Name[i], NPatlas$compound, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "NPatlas")
} else if (any(grepl(data.merged.filtered$GNPS_Compound_Name[i], mibig$compound, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "MiBIG")
} else if (any(grepl(data.merged.filtered$GNPS_smiles_nostereo[i], NPatlas$smiles_nostereo, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "NPatlas")
} else if (any(grepl(data.merged.filtered$GNPS_smiles_nostereo[i], mibig$smiles_nostereo, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "MiBIG")
} else if (any(grepl(data.merged.filtered$DEREP_Name[i], NPatlas$compound, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "NPatlas")
} else if (any(grepl(data.merged.filtered$DEREP_Name[i], mibig$compound, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "MiBIG")
} else if (any(grepl(data.merged.filtered$DEREP_smiles_nostereo[i], NPatlas$smiles_nostereo, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "NPatlas")
} else if (any(grepl(data.merged.filtered$DEREP_smiles_nostereo[i], mibig$smiles_nostereo, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "MiBIG")
} else if (any(grepl(data.merged.filtered$PCDL_Compound_Name[i], NPatlas$compound, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "NPatlas")
} else if (any(grepl(data.merged.filtered$PCDL_Compound_Name[i], mibig$compound, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "MiBIG")
} else if (any(grepl(data.merged.filtered$CSI_name[i], NPatlas$compound, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "NPatlas")
} else if (any(grepl(data.merged.filtered$CSI_name[i], mibig$compound, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "MiBIG")
} else if (any(grepl(data.merged.filtered$CSI_smiles[i], NPatlas$smiles_nostereo, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "NPatlas")
} else if (any(grepl(data.merged.filtered$CSI_smiles[i], mibig$smiles_nostereo, ignore.case = TRUE), na.rm = TRUE)) {
  df <- rbind(df, "MiBIG")
}  else {
  df <- rbind(df, NA)
}
i <- i+1
print(i)
print(data.merged.filtered$GNPS_Compound_Name[i])
}

data.m.f.annotated <- cbind(data.merged.filtered, df) #updates metadata table
colnames(data.m.f.annotated)[27] <- "npStatus"

##Propagate by gnps cluster index

gnpsClusters <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/Annotations/GNPS-clusters.csv")
data.clustered <- merge(gnpsClusters, data.m.f.annotated,
                        by = 'feature.ID', all = TRUE)

i <- 1
while (i <= nrow(data.clustered)) {
  if (!is.na(data.clustered$npStatus[i]) & data.clustered$componentindex[i] > -1) {
    CI <- data.clustered$componentindex[i]
    j <- 1
    while (j <= nrow(data.clustered)) {
      if (data.clustered$componentindex[j] == CI) {
        data.clustered$npStatus[j] <- "DB hit/cluster"
      }
      j <- j + 1
    }
  }
  i <- i + 1
}

##Append NPClassifer

NPCIDs <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/Annotations/NPClassifier.csv")
NPCIDs <- separate(NPCIDs, col = id, into = c("siriusID", "sirius", "feature.ID"), sep = "\\_")
NPC.form <- NPCIDs[, 3:4]
NPC.prob <- filter(NPCIDs, NPC.superclass.Probability >= 0.7 & NPC.superclass.Probability != "N/A")  #Accept prob. > 70%
NPC.prob <- NPC.prob[, c(3, 8:11)]
colnames(NPC.prob) <- paste("NPC", colnames(NPC.prob), sep = "_")
colnames(NPC.prob)[1] <- "feature.ID"
colnames(NPC.form) <- paste("NPC", colnames(NPC.form), sep = "_")
colnames(NPC.form)[1] <- "feature.ID"

data.clust.NPC <- merge(data.clustered, NPC.form,
                        by = 'feature.ID', all = TRUE)
data.clust.NPC <- merge(data.clust.NPC, NPC.prob,
                        by = 'feature.ID', all = TRUE)

##Filter out primary metabolites

data.clust.filtered <- filter(data.clust.NPC, is.na(GNPS_npclassifier_class) |
                              GNPS_npclassifier_class != "Unsaturated fatty acids" &
                             GNPS_npclassifier_class != "Sphingoid bases" &
                             GNPS_npclassifier_class != "Purine nucleos(t)ides" &
                               GNPS_npclassifier_class != "Saccharides" &
                               GNPS_npclassifier_class != "Monosaccharides" &
                               GNPS_npclassifier_class != "Dipeptides" &
                               GNPS_npclassifier_class != "Aminoacids|Dipeptides" &
                               GNPS_npclassifier_class != "Aminoacids" &
                               GNPS_npclassifier_class != "Dipeptides|Tripeptides" &
                               GNPS_Compound_Name != "Ursolic acid" &
                               GNPS_subclass != "Purines and purine derivatives")
data.clust.filtered <- filter(data.clust.filtered, is.na(NPC_NPC.superclass) |
                               NPC_NPC.superclass != "Amino acid glycosides" &
                               NPC_NPC.superclass != "Aminosugars and aminoglycosides" &
                               NPC_NPC.superclass != "Fatty Acids and Conjugates" &
                               NPC_NPC.superclass != "Fatty amides" &
                               NPC_NPC.superclass != "Fatty esters" &
                               NPC_NPC.superclass != "Glycerolipids" &
                               NPC_NPC.superclass != "Glycerophospholipids" &
                               NPC_NPC.superclass != "Saccharides" &
                                NPC_NPC.superclass != "Nucleosides" &
                                NPC_NPC.superclass != "Apocarotenoids" &
                                NPC_NPC.superclass != "Fatty acyl glycosides" &
                                NPC_NPC.superclass != "Fatty acyls" &
                                NPC_NPC.superclass != "Nicotinic acid alkaloids" &
                                NPC_NPC.superclass != "Eicosanoids" &
                                NPC_NPC.superclass != "Sphingolipids" &
                                NPC_NPC.superclass != "Flavonoids" &
                                NPC_NPC.superclass != "β-lactams")
data.clust.filtered <- filter(data.clust.filtered, is.na(NPC_NPC.class) |
                                NPC_NPC.class != "Aminoacids" &
                                NPC_NPC.class != "Aminosugars" &
                                NPC_NPC.class != "Glycerophosphocholines")

app.feats <- data.clust.filtered
app.feats <- select(app.feats, -row.ID)

featureIDs <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/feature_id_to_mass_pre-tidied.csv")

app.feats <- merge(featureIDs, app.feats, by = "feature.ID")
#app.feats is used in script_pca
