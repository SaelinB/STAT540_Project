# Changes

## Dataset
We are still using the same dataset, which contains paired bulk and single cell RNAseq.

## Changes in Analysis
After performing initial differential expression analysis between bulk and bulkized data across all samples, we found that around 83% genes are differentially expressed. This is unexpected and suggests some underlying issue, since we don't expect the data to be so significantly different. Our hypothesis is that there are significant differences in expression levels between samples between bulk/bulkized data such that it confounds this type of analysis. To tackle this issue, we've made two attempts. First, we tried to correct some expression level biases by regressing out differences. Second, we are attempting to identify "interesting" genes as regression outliers between bulk and bulkized per sample. This effectively subverts the between-sample issue. Afterwards, genes that are identified and occur in multiple samples (above some to-be-determined threshold number) will be followed up on.

## Changes in Assignments
In addition to the previously assigned tasks, Vladimir has attempted to correct the bulkized data for any biases in the protocol by fitting a regression line between bulk and bulkized per sample and adjusting the bulkized data so that the coefficient of this line is ~1 and the intercept is ~0. The motivation behind this is to equalize expression levels, on average, between bulk/bulkized for each sample so that the data is directly comparable.

Jordan has adopted an assignment in retrieving marker genes from Neuroexpresso and associating genes with data from BiomaRt. He has also attempted to see whether there is a bias in genes found to be differentially expressed towards specific cell types by comparing Neuroexpresso marker genes to the differentially expressed genes that are also marker genes to see whether there is any discrepancy.

# Analysis Progress
Much analysis has been done, but results are still far from final. Most importantly, the general workflows for the final analyses we want to do have been put in place and are being worked on to extract more desirable results.

## Data Preprocessing
Data has been effectively preprocessed as described in the original proposal. Marker genes, bulkized counts, BiomaRt data and plottable long counts are all saved for downstream analysis. Correlation between bulk and bulkized counts has been assessed.

## Data Quality/Cleaning
Data and metadata have both been inspected for consistency and quality. P-value distributions of differentially expressed genes have been inspected. Genes with low expression have been dropped in analysis and low quality cells (high fraction of mitochondrial RNA or ribosomal RNA) have been confirmed to be ignored.

## Differential Expression Analysis
Some preliminary differential expression analysis has been done. An analysis pipeline is in place, but will have to be adapted to deal with the problems outlined above.

## Marker Gene Analysis
Preliminary analysis of differentially expressed genes has been done to identify cell type, GC content, length and chromosome biases. However, these analyses are still rather rudimentary and will have to be improved to yield more statistically grounded results.

# Results
Preliminary results show more differentially expressed genes than we expect (approximately 30,000 genes). Of these differentially expressed genes, inspection doesn't reveal any obvious cell type, GC content, length or chromosome biases.

Correlation between expression profiles within individuals (bulk vs. bulkized) has been assessed and is very strong (r^2 > 0.8), indicating that in general, bulkized single cell RNAseq UMI counts and bulk RNAseq counts are quite comparable.
