## Use case
For converting NCBI datasets downloaded from [GDSbrowser](https://www.ncbi.nlm.nih.gov/sites/GDSbrowser) to CSV files. Here is a list of tasks::

 - Clean the data by removing all extra information
 - Extract genes
 - Extract chromosomes
 - Merge the repetitive genes by calculating average
 - Remove fully NULL columns
 - Remove columns with ####_at_ names
 - Remove columns with --Control names
 - Impute the data to fill the NAs
 - Normalize samples of each class based on their median
 
This code extracts data with and without features' names, features' names, chromosomes names and store them into four CSV files, respectively.

## To run

### Step 0
Download the *DataSet full SOFT file* from the [GDSbrowser](https://www.ncbi.nlm.nih.gov/sites/GDSbrowser) and extract it in a new folder

### Step 1
Download and copy this code into the same folder

### Step 2
Run the code in RStudio or R environment

### Step 3
Four new CSV files will be created in the same folder containing data with and without features' names, features' names, and chromosomes names
