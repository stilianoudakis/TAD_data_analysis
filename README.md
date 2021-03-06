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
   + 06.1mourad_model.Rmd - code for running the multiple logistic regression model proposed by Mourad et al.
   + 07_comparing_models - code for evaluating the performance of mourad model vs random forest and gbm using data ran through pipeline
   + `data_exploration` - folder with files that explore the distribution of each feature in the data
   + `evaluation_normalization` - folder with files exploring how normalizing the data (with and without log2 transformation) affects each model
   + `evaluating_SMOTE` - folder with files exploring the SMOTE function from the `DMwR` package, including which combination of perc.over and perc.under performs best with the data
   + `evaluating_variable_reduction` - 
      + `Boruta` - using the boruta package to reduce the features of the data (applied to both the full data and a balanced dataset from randomly sampling and comparing)
      + `variable_selection` - using stepwise selection to filter the feature space (forward, backward, and both) on a balanced dataset using both random sampling and SMOTE and comparing
      + `RFE` - variable selection using recursive feature elimination

* `K562`  
   + 01arrowhead_data
   + 02TAD_boundary_analysis
   + 03incorporating_distance
   + 04adding_annotations
   + 05model_filtering
   + variable_selection.R
   + 06.1mourad_model.Rmd 
   + 07model_building
   + 08model_performance
   + `data_exploration`

* `reports`  
   + chr1_report - report evaluating elastic net, random forest, gbm, and glm models on chr1 TAD boundary data
   + data_exploration - univariate analyses on features for GM12878 cell line
   + data_explotation_k562 - univariate analyses on features for k562 cell line
   + investigating_tad_boundaries - 
   + normalization_report - evaluating different normalization techniques
   + smote_report - evaluating different class balancing techniques
   + variable_selection_report - evaluating different variable selection techniques
   + measuring_performance_normalization - tables and figures measuring performance of different normalization techniques
   + measuring_performance_reduction - tables and figures measuring performance of different variable selection techniques
   + measuring_performance_smote - tables and figures measuring performance of different class balancing techniques
   + comparing models - comparing performances of the mourad models (with and without LASSO estimation) with random forest using data ran through pipeline
   
  
* `presentations`
   + SSTP_2018.Rmd

* Papers - Important and relevant papers to read

* Project - Overview and objectives of the current project

* ToDo - list of objectives to complete accompanied by dates

* misc - Miscellaneous code saved for later use


