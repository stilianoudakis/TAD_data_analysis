# TAD_data_analysis
Analysis of Topologically Associated Domain data and their respective boundaries


## Folders & Files

* `data`  
   + arrowhead_data.txt - The TAD boundary data that was retrieved from from the paper by Rao, Huntley, et. al titled, "A three-dimensional map of the human genome at kilobase resolution reveals prinicples of chromatin looping" for the GM12878 cell line (GEO accession: GSE63525). 
   + arrowhead_k562_data.txt - Same data as `arrowhead_data.txt` except for the K562 cell line.
   + GSE63525_GM12878_subcompartments.bed - subcompartment data classifying genomic coordinates as being in specific subcompartments A and B.

* `GM12878` 
   + 01arrowhead_data - cleaning and saving the TAD boundary data for downstream analysis
   + 02TAD_boundary_analysis - creating 1kb TAD boundary bins and creating variables for whether or not the bin overlaped a particular subcompartment
   + 03incorporating_subcomp_distance - adding distance variables for each of the two subcompartment variables (A and B)
   + 04adding_annotations - creating features from genomic annotations in the form of a binary variables and continuous distance variables for each annotation
   + 05model_filtering - filtering the data by removing features with near zero variances and turning binary variables to factors
   + 06model_building - code for various maching learning algorithms
   + 07model_performance - code for evaluating the performance of each algorithm
   + `data_exploration` - folder with files that explore the distribution of each feature in the data
   + `evaluation_normalization` - folder with files exploring how normalizing the data (with and without log2 transformation) affects each model
   + `evaluating_SMOTE` - folder with files exploring the SMOTE function from the `DMwR` package, including which combination of perc.over and perc.under performs best with the data
   + `evaluating_variable_reduction` - 
      + `Boruta` - using the boruta package to reduce the features of the data (applied to both the full data and a balanced dataset from randomly sampling and comparing)
      + `variable_selection` - using stepwise selection to filter the feature space (forward, backward, and both) on a balanced dataset using both random sampling and SMOTE and comparing

* `K562`  
   + 01arrowhead_data
   + 02TAD_boundary_analysis
   + 03incorporating_distance
   + 04adding_annotations
   + 05model_filtering

* `reports`  
   + chr1_report
   + data_exploration
   + investigating_tad_boundaries
   + normalization_report
   + smote_report

* Papers - Important and relevant papers to read

* Project - Overview and objectives of the current project

* ToDo - list of objectives to complete accompanied by dates

* misc - Miscellaneous code saved for later use


