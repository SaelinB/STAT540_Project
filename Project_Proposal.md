# Motivation
The advent of single-cell RNAseq (scRNAseq) has allowed for the understanding of gene expression at an unprecedented level. As such, it has been an area of high interest. However there remain many technically difficult hurdles to overcome, not the least of which being sample preparation.

Conventional single-cell transcriptomic preparation protocols usually involve the following steps [[1]](#1):
1. Dissection of whole tissue
2. Enzymatic breakdown of ECM
3. Optional cell type-specific enrichment

Careful handling and processing of cells destined for scRNAseq is important for preserving the native expression profile. Due to the sensitivity of these cells, recent advances in bioengineering has led to the development of some microfluidic devices that are optimized for performing these steps. However, it is still common for cell preparation to be done by hand [[2]](#2).

When cells are subjected to undue stress, rupture or death, it is common to term them "low quality" [[3]](#3) as they are mainly unsuitable for analysis. While some such low quality cells never make it to the sequencer in scRNAseq experiments, those that do are typically filtered out in downstream analysis by a variety of methods [[4](#4),[5](#5)]. This raises an important question: To what degree is the single-cell RNA profile representative of the whole tissue? If cell type-specific enrichment, (ie. via flow cytometry or differential centrifugation) is not performed, it is often taken for granted that the underlying cellular population is unchanged by the preparation protocol [[6]](#6). In actuality, this may not be the case since, for example, some cell types (ie. neurons) are more easily disrupted than others by standard tissue dissociation protocols [[7]](#7) and may skew data towards more "durable" cell types. Thus, it is not unreasonable to believe that large-scale preparation-induced cell quality issues do not arise stochastically. To illustrate this, recent studies examining the correlation between cell type proportions (arising from either from scRNAseq or _in situ_ counts) and estimates from bulk tissue (using marker gene signatures) reveal, at best, a Pearson correlation coefficient of only approximately 0.4 [[8]](#8) (a moderately weak correlation) and some significant differences in estimated absolute cell type proportions between methods [[8]](#8) indicating that either (a) these proportions are only moderately predictable using marker gene signatures or (b) the cellular makeup is, in fact, different between the two. If the latter situation is the case, elucidating such biases will be useful in ensuring scRNAseq preparation protocols are amenable to the tissue under study. By nature of unavoidable differences in exact sequencing protocols (ie. for different platforms, by different personnel, etc.), the results of our analysis may not offer specific cell type proportion "correction factors" for various other scRNAseq studies, but it nonetheless takes a step towards uncovering yet unknown limitations in scRNAseq interpretability.

To study these hypothesized technical biases, we will make use of recently published datasets that provide paired bulk and single cell RNAseq from snap-frozen deceased human brain samples (outlined below). These datasets draw from paired sample extractions such that each sequenced tissue sample is derived from the same tissue sample at the same developmental time and spatial region. As such, it is a fair assumption that each replicate extraction is nearly identical at the time of harvesting. Thus, if we detect genes that exhibit differential expression (DE) between these two datasets, we can attribute this signal, in large part, to technical factors (ie. any event that happened from the time of harvesting until quantification, which is necessarily non-biological). It is important here to note that some bias in quantification may also arise from sequencing-derived bias, which we can not rule out as an additional confound. Once genes exhibiting some DE between the scRNAseq and bulk tissue RNAseq datasets are found, we will analyze them for enrichment of specific classical cell types and other data-driven groupings (ie. functional enrichment analysis).

# Previous Work
While the idea of cell quality in scRNAseq has been heavily discussed, this specific idea has not yet been explored in published literature that we know of.

# Division of Labour
| Name | Department/Program | Expertise/Interests | GitHub ID | Responsibilities |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| Lucía Quesada-Ramírez | GSAT | Plant genomics, polyploidy | @luciaq21 | Poster design, Differential expression analysis |
| Saelin Bjornson | Bioinformatics | Comparative genomics, phylogenetics, evolution | @SaelinB | Enrichment analysis, PCA/RRHO and metadata quality checking |
| Vladimir Nikolic | Bioinformatics | Genome assembly, algorithms | @schutzekatze | Enrichment analysis, data validation, data visualization |
| Jordan Sicherman | Bioinformatics | snRNAseq, differential expression analysis | @jsicherman | Data processing/preparation, project planning, data quality checking |

# Datasets

We have access to the following datasets which would be of use in our analysis: [Single-cell genomics identifies cell type-specific molecular changes in autism](https://science.sciencemag.org/content/364/6441/685).

The single cells (approximately 120,000) in this dataset originate from 48 post-mortem tissue samples from the prefrontal, anterior cingulate and insulate cortical regions of 32 patients (16 with ASD, 16 control). They were sequenced on the 10x Genomics Chromium platform and demultiplexed, aligned and quantified using 10x Genomics CellRanger software [[9]](#9). This data was further filtered and quality checked, as outlined in the [data description](data/raw/README.md). This data also contains paired RNAseq data of bulk tissue taken from tissue adjacent to that used for single cell isolation.

We additionally plan to use a published set of brain cell type marker genes from a meta-analysis done in [NeuroExpresso](https://www.eneuro.org/content/4/6/ENEURO.0212-17.2017) for performing cell type enrichment. As these are mouse-derived marker genes, we will use [HomoloGene](https://www.ncbi.nlm.nih.gov/homologene) to detect their human homologs.

# Specific Aims and Methodology
1. "Bulkize" the single cells from our scRNAseq data.
  - "Bulkizing" refers to the straightforward sum of gene expression among all cells of a given sample.
  - After bulkization, these counts should correlate with counts in the real bulk RNAseq samples. This will be assessed using the Pearson correlation coefficient.
2. Perform DE analysis between bulkized single cells and their matched bulk tissue data.
  - Typically, Limma is inappropriate for scRNAseq DE analysis (since, in short, the single cells don't follow the expected distribution). However, "bulkizing" the cells should resolve this issue - this will also be examined by inspecting the distribution of counts.
 - Genome-wide DE analysis will be done using Limma-voom using a standard linear model. As each of the 48 comparisons will be done on the same sample (with the same metadata), no other covariates besides the categorical groups "bulk" or "bulkized" are needed. 
3. Perform enrichment analysis on our genome-wide DE gene hits from (2) using a published set of cell-type specific marker genes.
  - PCA will be done on our gene profiles, to which we will map our marker genes. This will confirm that cell types can be inferred by our marker gene set. 
  - Genes of our marker gene set will be extracted and ranked by strength of DE. 
  - Standard rank-rank hypergeometric overlap of our marker gene set between bulk and bulkized groups will identify the strength and pattern of correlation between the two groups. In addition, standard rank-rank hypergeometric overlap can be done between the differentially-expressed gene set and the marker gene set to determine if differentially expressed genes are enriched in marker genes. 
4. Time permitting, perform GO Enrichment on DEG hits from (2) to see if our gene set has an enrichment of additional functional annotations not included in our marker gene set (for example, as a result of expression differences due to different protocols). 
  - We will use Broad Institute's GSEA tool to identify representative pathways with positive correlations in our gene set.
  - These pathways will be compared to literature to further understand and validate our findings.
5. Time permitting, to validate our findings, we will attempt to correct for any cell type biases we detect in the scRNAseq data and perform sample-level correlation of corrected cellular proportions in the scRNAseq data and estimated cellular proportions from the paired bulk tissue RNAseq.
  - Cell type biases will be corrected ad hoc, based on the magnitude of the expression differences.
  - Correlation will be assessed as previously explored in literature [[8]](#8).

# References
<a id="1">[1]</a> 
Nguyen, Q.H., Pervolarakis, N. and Kessenbrock, K. (2018). 
Experimental Considerations for Single-Cell RNA Sequencing Approaches. 
Frontiers in Cell and Developmental Biology, 6(108).

<a id="2">[2]</a> 
Qiu, X., Jesus, J.D., Pennell, M., Troiani, M. and Haun, J.B. (2014). 
Microfluidic device for mechanical dissociation of cancer cell aggregates into single cells. 
Lab on a Chip, 15(1).

<a id="3">[3]</a> 
Ilicic, T., Kim, J.K., Kolodziejczyk, A.A. et al. (2016). 
Classification of low quality cells from single-cell RNA-seq data. 
Genome Biology, 17(29).

<a id="4">[4]</a> 
Brennecke, P., Anders, S., Kim, J.K. et al. (2013). 
Accounting for technical noise in single-cell RNA-seq experiments. 
Nature Methods, 10, 1093-1095.

<a id="5">[5]</a> 
Islam, S., Zeisel, A., Joost, S. et al. (2014). 
Quantitative single-cell RNA-seq with unique molecular identifiers. 
Nature Methods, 11, 163-166.

<a id="6">[6]</a> 
Mancarci, B.O, Toker, L., Tripathy, S.J., Li, B., Rocco, B., Sibille, E. and Pavlidis, P. (2017). 
Cross-Laboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk Tissue Data. 
eNeuro, 4(6).

<a id="7">[7]</a> 
Habib, N., Li, Q., Regev, A. et al. (2016). 
Div-Seq: Single-nucleus RNA-Seq reveals dynamics of rare adult newborn neurons. 
Science, 353(6302), 925-928.

<a id="8">[8]</a> 
Patrick, E., Taga, M., Mostafavi, S. et al. (2019). 
Deconvolving the contributions of cell-type heterogeneity on cortical gene expression. 
bioRxiv.

<a id="9">[9]</a> 
Velmeshev, D., Schirmer, L., Jung, D. et al. (2019). 
Single-cell genomics identifies cell type-specific molecular changes in autism. 
Science, 364(6441), 685-689.
