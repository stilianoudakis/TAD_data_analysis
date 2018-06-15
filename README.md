# TAD_data_analysis
Analysis of Topologically Associated Domain data and their respective boundaries

- `Papers.md` - paper notes

- `Project.md` - project draft

## data

* arrowhead_data(2).txt: The TAD boundary data used for analyses cleaned from the GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt file. (2) refers to the same data with quotations removed around the Chromosome column

* subcompartments.bed: subcompartment data classifying genomic coordinates as being in specific subcompartments A and B.

## `01arrowhead_data.Rmd`

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

## `02investigating_tad_boundaries.Rmd`

Investigating the different type of boundaries that occur in the domain data.

* Purpose: To investigate the different ways TAD boundaries share coordinates
* The Data: arrowhead_data.txt (GSE63525_GM12878_primary+replicate_Arrowhead_domainlist)
* Methods: Various data manipulation techniques
* Take Home Message: The are 4 possible was TAD boundaries can share a coordinate:
   1. Same Start site
   2. Same End site
   3. Consecutive End and Start sites
   4. Duplicate intervals (none of these)

## `03TAD_boundary_analysis.Rmd`

* Purpose: To create flanked TAD boundary regions of size 1kb (500 bases on each side of TAD boundary) and identify where on the 1 kb binned genome they appear. Likewise, identifying what compartment each of these bins correspond to.
* The Data: arrowhead_data.txt, GSE63525_GM12878_subcompartments.BED
* Methods: 
   1. A vector of TAD boundaries was created by concatenating the start and end coordinates of the arrowhead data.
   2. The vector was sorted and duplicates corresponding to shared coordinates were removed.
   3. The genome was binned at 1kb starting from the (minimum coordinate - 500) to the (maximum coordinate + 500) for each chromosome. Therefore, a TAD boundary will either appear in the middle of a bin or it will not.
   4. A vector Y was created that identified whether or not a TAD boundary was located in the bin or not
   5. Genomic annotation data was used in the form of subcompartments (A or B). A vector X was created which labled which subcompartment each 1kb bin fell into (A, B, or N for none). 
   6. A data frame containing the class Y and predictor X was created for the logistic regression model.
* Notes: Only overlaps fully contained in the subcompartment intervals were used. Therefore, some bins that partially overlapped a subcompartment will be labled N. Going further, these need to be labeled with the percentage of there overlap. Also, chromosome X was associated with no subcompartment (all bins were labeled N).


## chr1_logistic_regression.R

R code for running a logistic regression model using subcompartment coordinates as predictors


