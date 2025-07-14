Code used to generate results, figures and tables for Science submission: *Human land use drove persistence of Holocene plant assemblages*

The data used to generate this publication were obtained from the Neotoma Paleoecological Database: details are located in *1_dataAcquisitionProcessing.R*.

Running scripts will require a copy of the entire repository, including subfolders, which contain raw data and supplementary functions. Please add correct working directory entries in toplevel .R scripts and relative file paths will do the rest.

R scripts make extensive use of code-folding functionality in RStudio.

REPOSITORY TREE:

*./1_dataAcquisitionProcessing.R*: commented code to access Neotoma data and harmonize taxonomy. Raw outputs have been included in repository so this script does not need to be re-run to reproduce results.

*./2_chronologyModelling.R*: commented code to run rBacon chronology models on dates stored in Neotoma. This script also adjusts pollen counts using relative-pollen-production estimates. Raw outputs have been included in repository so this script does not need to be re-run to reproduce results. Note that due to size restrictions, rBacon visualizations and raw output files have been removed from the *./baconRuns/* subfolder.

*./3_ratepolEstimation.R*: commented code to run RRatePol paleoecological resampling. Modified RRatePol functions that allow this script to produce dissimilarities used for estimating ecological persistence are included in the *./functions/* subfolder. This includes Late Glacial and Holocene temporal extents and grains, as well as variant family-level analyses.

*./4_analyses.R*: commented code to produce all main text and supplementary figures and results. This includes Late Glacial and Holocene temporal extents and grains, as well as variant family-level analyses.

**./baconRuns**: Folder to store rBacon chronology model outputs. Each subfolder is a Neotoma datasetid, containing the specific chronology model. Note that due to size restrictions, rBacon visualizations and raw output files have been removed. If you intend to run *./2_chronologyModelling.R*, you will need to delete the entries in this folder.

**./functions**: Small ease-of-use functions read into top-levels scripts and modified RRatePol functions.

**./outputs**: File outputs from R scripts stored as ".rds" compressed file types. Includes intermediate outputs produced by each script so that *./4_analyses.R* can be run to reproduce results. Files with "trajectory" are outputs from *./3_ratepolEstimation.R*. Files with "modelData" or "holoData" are final raw data products from *./4_analyses.R* after subsetting and PCA processing.

**./plots**: Figure outputs from *./4_analyses.R*, used in main text and online supplement. Variants labelled "EDIT" have had labels moved in post-R software, including .png edits to include in submitted manuscript.

**./rawdata**: Input data objects and processed files from *./1_dataAcquisitionProcessing.R*. 