# Black Apple Project
The repository here contains code used in the data analysis and preparation of the manuscript: XXXX

## Tidying of annotations (Tidy-annotations.R)

The provided script is a data processing pipeline that involves various data manipulation and tidying operations to merge and filter data from multiple CSV files. Below is a step-by-step explanation of the script:

1. Loading of libraries: including dplyr, tidyr, stringr, and jsonlite
2. Importing, extraction and tidying of NP Databases: NPatlas and MIBIG
3. Merging annotations from GNPS, dereplicator, CSI-FingerID and PCDL in-house database with the appropriate metadata
4. Tidying compound names and removal of stereochemistry from SMILES in the annotation list and NP databases
5. Querying annotation list against the NP databases
6. Propagation of natural product classes according to GNPS clusters based on IDs in NP databases
7. Appending NPClassifier predictions
8. Filtering out primary metabolites
9. The script results in a dataframe "app.feats" which contains all the features across the black apples, and OSMAC studies of fungal isolates from black apples.

## Data analysis (Data-analysis.R)

Rather than working on the annotations, this script uses those annotations and the raw data output following MZMine to perform some data analysis and visualization using the dplyr, tidyr, and ggplot2 libraries in R. Below is a step-by-step explanation of the script:

1. Loading of libraries" dplyr, tidyr, viridis, and ggplot2
2. Importing of tidied GNPS quant file (GNPS_tidied.csv) and metadata (same format as GNPS is fine)
3. Removing primary metabolites using the tidy annotations (app.feats)
4. Tidying and merging of data with metadata
5. Removal of features associated with the blanks based on certain thresholds, creating the dataframe "data.filt.blanks"
6. Comparison of features associated with whole black apple extracts (Apple_natural), excisions of fungi from black apples (Apple_excision) or  fungal isolates cultivated in an OSMAC approach (Laboratory)
7. Extraction of features found in both the fungal isolates grown in the lab and on black apples (in.both2)
8. Visualisation of compound classes of fungal features found in black apples using NPclassifier superclasses. 

## Bubble chart (Bubble-chart.R)

This script generates a bubble chart visualising the proportion of different superclasses of fungal metabolites present in various black apple samples. Below is a step-by-step explanation of the script:

1. Loading of libraries: ggplot2 and reshape2
2. In a similar manner to that used for step 7 in data analysis above, the features associated with black apple fungi grown in the lab and features associated with black apples were extracted
3. Calculates the NP superclasses associated with each indvidual black apple extract
4. Bubble chart is created using the generated dataframe extracted



Calculating Superclasses Proportion: The code groups the features detected in both fungus and apple samples based on their superclasses. It calculates the proportion of each superclass and prepares the data for the bubble chart.

Generating Bubble Chart: The code iterates through different apple samples and calculates the proportion of superclasses for each sample using the apple.superclasses function. It then combines the data into a new dataframe df3.

Merging with Apple IDs: The code reads another CSV file, "appleIDs.csv," containing information about apple IDs, and merges it with the dataframe df3 to convert filenames to appleIDs.

Plotting the Bubble Chart: Finally, a bubble chart is created using ggplot2, visualizing the proportion of different superclasses for each apple sample. The size of each bubble represents the proportion of metabolites in that superclass, and the color of the bubbles represents the superclass category.
