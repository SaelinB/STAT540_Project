library(homologene)
library(biomaRt)

# Get human homologs ---------------
mus_genes <- read.csv('data/raw/neuroexpresso.cortex.tsv', sep = '\t')
hum_genes <- apply(mus_genes, 2, function(x) mouse2human(x)$humanGene)

# Convert to Ensembl symbols ---------------
ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

# Return the HGNC symbol, Ensembl ID, GO ID, chromosome name, transcript length and GC %
getProperties <- function(geneList, filter) {
  getBM(mart = ensembl, attributes = c('hgnc_symbol',
                                       'ensembl_gene_id',
                                       'go_id',
                                       'chromosome_name',
                                       'transcript_length',
                                       'percentage_gene_gc_content'), filters = filter, values = geneList)
}

# Do it for marker genes
marker.properties <- lapply(hum_genes, getProperties, filter = 'hgnc_symbol')

# Do it for all the genes of interest
all.genes <- unique(readRDS('data/processed/counts.merged.long.rds')$gene)
all.properties <- getProperties(all.genes, 'ensembl_gene_id')

# Save ---------------
saveRDS(marker.properties, 'data/processed/genes.markers.rds')
saveRDS(all.properties, 'data/processed/genes.detected.rds')