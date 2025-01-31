# Black Apple Project
The repository here contains code used in the data analysis and preparation of the manuscript: Unveiling the Fungal Diversity and Associated Secondary Metabolism on Black Apples ([https://doi.org/10.1101/2023.11.02.565319](https://doi.org/10.1128/aem.00342-24)). The GNPS-FBMN workflow and results can be found here: ([https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=71057036992e4ac6aef647aab6cf7557](https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=71057036992e4ac6aef647aab6cf7557)), and the associated raw data under MASSIVE MSV000092823.

If you use or ammend any of the scripts utilised in this work, please cite the above publication.

## Tidying of annotations (Tidy-annotations.R)

The provided script is a data processing pipeline that involves various data manipulation and tidying operations to merge and filter data from multiple CSV files. Below is a step-by-step explanation of the script:

1. Loading of libraries: including dplyr, tidyr, stringr, and jsonlite
2. Importing, extraction and tidying of NP Databases: NPAtlas and MiBIG
3. Merging annotations from GNPS, dereplicator, CSI-FingerID and PCDL in-house database with the appropriate metadata
4. Tidying compound names and removal of stereochemistry from SMILES in the annotation list and NP databases
5. Querying annotation list against the NP databases
6. Propagation of natural product classes according to GNPS clusters based on IDs in NP databases
7. Appending NPClassifier predictions
8. Filtering out primary metabolites
9. The script results in a dataframe "app.feats" which contains all the features across the black apples, and OSMAC studies of fungal isolates from black apples.

## Data analysis (Data-analysis.R)

Rather than working on the annotations, this script uses those annotations and the raw data output following MZMine to perform some data analysis and visualization using the dplyr, tidyr, and ggplot2 libraries in R. Below is a step-by-step explanation of the script:

1. Loading of libraries: dplyr, tidyr, viridis, and ggplot2
2. Importing of tidied GNPS quant file (GNPS_tidied.csv) and metadata (same format as GNPS is fine)
3. Removing primary metabolites using the tidy annotations (app.feats)
4. Tidying and merging of data with metadata
5. Removal of features associated with the blanks based on certain thresholds, creating the dataframe "data.filt.blanks"
6. Comparison of features associated with whole black apple extracts (Apple_natural), excisions of fungi from black apples (Apple_excision) or  fungal isolates cultivated in an OSMAC approach (Laboratory)
7. Extraction of features found in both the fungal isolates grown in the lab and on black apples (in.both2)
8. Visualisation of compound classes of fungal features found in black apples using NPclassifier superclasses.

## Superclass piecharts (Superclass-pie-chart.R)

This script goes a level deeper than that for step 8 in Data-analysis.R, and can profile the superclasses for a given fungal genera (with the potential to generate based on species, and/or other criteria).

## Bubble chart (Bubble-chart.R)

This script generates a bubble chart visualising the proportion of different superclasses of fungal metabolites present in various black apple samples. Below is a step-by-step explanation of the script:

1. Loading of libraries: ggplot2 and reshape2
2. In a similar manner to that used for step 7 in data analysis above, the features associated with black apple fungi grown in the lab and features associated with black apples were extracted
3. Calculates the NP superclasses associated with each indvidual black apple extract using the created apple.superclasses
4. Black apple extract IDs are appended with an external metadata file "appleIDs.csv"
5. Bubble chart is created using ggplot2, visualizing the proportion of different superclasses for each apple sample. The size of each bubble represents the proportion of metabolites in that superclass, and the color of the bubbles represents the superclass category

## Fungal species present in given black apple extracts (Species-checker.R)

This script generates a heatmap and a barchart indicating the frequency of species present in given black apple extracts. Below is a step-by-step explanation of the script:

1. Extracts a list of features indicative of a specified species (e.g. "polonicum")
2. Queries black apple extracts using the apple.finder functional to find which contain features associated with the above list of features (with strict qualification criteria)
3. Using the above code, a user-made .csv file was made to tabulate which black apple extracts contained which fungal species
4. Generates a heatmap using this user-made .csv file indicating the presence of each species in each individual black apple extract
5. Generates a barchart indicating the frequencing each species is detected in the black apples as a whole

## Features present in given black apple extracts (Feature-checker.R)

This script can be used to query the frequency of a given feature or features present in given black apple extracts. Below is a step-by-step explanation of the script:

1. User input features to query
2. Sets up the functions needed to query black apples
3. Generates a data frame of which apples contain the queried features
4. Calculates the number of black apples which contain at least one of the queried features

## Features present in given black apple extracts (Perc-Abundances.R)

This script calculates the average abundance of features present in black apples compared to the OSMAC study.  Below is a step-by-step explanation of the script:

1. Normalises data in each sample to the maximum peak area
2. Calculates the average abundance for each feature across the apple extracts
3. Calculates the average abundance for each feature across the OSMAC extracts
4. Merges the data together
5. Calculates the ratio of the apple averages to the OSMAC averages
