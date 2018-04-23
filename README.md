# TAD_data_analysis
Analysis of Topologically Associated Domain data and their respective boundaries

## arrowhead_data.Rmd

Cleaning and changing the characteristics of the TAD boundary data.

The contact domain annoation file is named: GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt

The arrowhead algorithm was used to determine domains

The columns are defined as follows: 
chromosome1, x1, x2, chromosome2, y1, y2, color, corner_score, Uvar, Lvar, Usign, Lsign

Explanations of the fields are as follows:

* chromosome = the chromosome that the domain is located on

* x1,x2/y1,y2 = the interval spanned by the domain (contact domains manifest as squares on the diagonal of a Hi-C matrix and as such: x1=y1, x2=y2)

* color = the color that the feature will be rendered as if loaded in Juicebox 

## investigating_tad_boundaries.Rmd

Investigating the different type of boundaries that occur in the domain data.

## TAD_boundary_analysis.rmd

Analysis of TAD boundary data retrieved from from the paper by Rao, Huntley, et. al titled,
"A three-dimensional map of the human genome at kilobase resolution reveals prinicples of chromatin looping"

GEO accession: GSE63525

## chr1_logistic_regression.R

R code for running a logistic regression model using subcompartment coordinates as predictors


