# TAD_data_analysis
Analysis of Topologically Associated Domain data and their respective boundaries

## data

* arrowhead_data(2).txt: The TAD boundary data used for analyses cleaned from the GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt file. (2) refers to the same data with quotations removed around the Chromosome column

* subcompartments.bed: subcompartment data classifying genomic coordinates as being in specific subcompartments A and B.

## arrowhead_data.Rmd

The TAD boundary data was retrieved from from the paper by Rao, Huntley, et. al titled,
"A three-dimensional map of the human genome at kilobase resolution reveals prinicples of chromatin looping"

GEO accession: GSE63525

The contact domain annoation file is named: GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt

The columns are defined as follows: 
chromosome1, x1, x2, chromosome2, y1, y2, color, corner_score, Uvar, Lvar, Usign, Lsign

Explanations of the fields are as follows:

* chromosome = the chromosome that the domain is located on

* x1,x2/y1,y2 = the interval spanned by the domain (contact domains manifest as squares on the diagonal of a Hi-C matrix and as such: x1=y1, x2=y2)

* color = the color that the feature will be rendered as if loaded in Juicebox 

## investigating_tad_boundaries.Rmd

Investigating the different type of boundaries that occur in the domain data.

* Purpose: To investigate the different ways TAD boundaries share coordinates
* The Data: arrowhead_data.txt (GSE63525_GM12878_primary+replicate_Arrowhead_domainlist)
* Methods: Various data manipulation techniques
* Take Home Message: The are 4 possible was TAD boundaries can share a coordinate:
   1. Same Start site
   2. Same End site
   3. Consecutive End and Start sites
   4. Duplicate intervals (none of these)

## TAD_boundary_analysis.rmd

* Purpose: To create genomic bins and flank the TAD boundaries
* The Data: arrowhead_data.txt, GSE63525_GM12878_subcompartments.BED
* Methods: Initial focus was paid to chromosome 22
   1. The starting and ending genomic coordinates were determined from the 1kb contact matrix. They were 16,049,000 and 51,244,000
   2. Genomic bins of size 50 bases were created for chromosome 22
   3. For the TAD boundary data, looking at chromosome 22, the start and end sites were concatenated together, duplicates were removed. This results in a vector of TAD boundary sites.
   4. Each boundary site was flanked on either side by 500 bases for a 1kb centered boundary region
   5. A vector Y indicating whether the bins overlaped the boundary intervals was created
   6. Genomic feature data was imported in the form of subcompartment classifications (A or B). Predictor vectors corresponding to if there was an overlap for each compartment were created
* Notes: There are only 36 subcompartments associated with chromosome 22. With a mean width of 1366668 for subcompartment A and 614287 for subcompartment B. Therefore the intervals are much larger than the genomic bins of 50 bases. It appears all partial overlaps land on the very end coordinate of the genomic bins thus the percentage overlap would be 1/50=.02. 


## chr1_logistic_regression.R

R code for running a logistic regression model using subcompartment coordinates as predictors


