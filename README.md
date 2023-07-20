# Black Apple Project
The repository here contains code used in the data analysis and preparation of the manuscript: XXXX

## Tidying of annotations
[Tidy-annotations.R]([URL](https://github.com/michael-cowled/black-apples/blob/main/Tidy-annotations.R))

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

## Data analysis
SCRIPT.YYY

