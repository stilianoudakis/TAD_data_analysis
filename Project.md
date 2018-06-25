# Epigenomic signatures of topologically associated domains (TADs)

TADs are genomic regions with highly interacting regions within a TAD and much less interactions with other TADs. TAD boundaries are defined as genomic regions where one TAD ends and another begins. These TAD boundaries have been shown to be enriched in various histone marks. Our goal is to better understand epigenomic signatures of TADs and TAD boundaries.


## Input data

- Most of the data will be in BED format, genomic coordinates. Note that Hi-C data resolution is relatively coarse, from 1kb to 1Mb genomic regions
- TAD boundaries
- TAD centers. Justification - [@Boulos:2013aa] showed that TAD boundaries coincide with the highest mean replication timing (MRT). Consequently, centers have the lowest MRT, thus important.
- Epigenomic data (transcription factors, histone marks, DNAse hypersensitive sites, chromatin states), cell type specificity is important


## Methods

- Classification problem - distinguishing TAD boundaries from non-TAD regions. TAD boundaries from TAD centers. TAD centers from non-center regions
- Need a strategy to assign epigenomic annotations to Hi-C genomic bins
    - Binary - if an annotation overlaps a region, use indicator "1"
    - Counts - how many peaks are in a region?
    - Signal strength - e.g., partitioning epigenomic signal into 20 signal strength windows, as in [@Di-Pierro:2017aa]
- Machine learning methods, with cross-validation


## Open questions

- Can we factor in distance
- Need to investigate class disbalance problem




