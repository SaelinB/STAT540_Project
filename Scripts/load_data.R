library(dplyr)
library(Matrix)
library(reshape2)
library(tibble)
library(parallel)
options('mc.cores' = 16) # Use 16 cores

# Un-normalize ---------------
# Load the normalized count matrix and associated meta...
counts.norm <- readMM('data/raw/counts.singlecell.normalized.mtx') %>% as('dgCMatrix')
metadata <- read.csv('data/raw/meta.singlecell.csv')

# Back-transform the normalized count matrix to a proper count matrix.
# Mutiplication occurs column-wise and we want row-wise (across cells), so doubly transpose.
# At least on my machine, this doesn't fit in memory, so divide into chunks
counts.unlogged <- do.call(cbind, mclapply(cumsum(c(rep(6534, 15), 6549)), function(chunk) {
  2^counts.norm[, (chunk - ifelse(chunk == 104559, 6548, 6533)):chunk] - 1
})) %>% as.matrix %>% Matrix(sparse = T) # This stypidity bypasses some weird dense issue

counts.raw <- t(t(counts.unlogged) * metadata$UMIs) / 10000
rm(counts.unlogged)

# Bulkize ---------------
# For each sample, spit back row means and then merge the columns together.
counts.bulkized.raw <- do.call(cbind, mclapply(unique(metadata$sample), function(sample) {
  counts.raw[, which(metadata$sample == sample)] %>% Matrix::rowSums(na.rm = T)
})) %>% Matrix(sparse = T)

counts.bulkized.norm <- do.call(cbind, mclapply(unique(metadata$sample), function(sample) {
  counts.norm[, which(metadata$sample == sample)] %>% Matrix::rowSums(na.rm = T)
})) %>% Matrix(sparse = T)
rm(counts.raw, counts.norm)

# Save our files.
writeMM(counts.bulkized.raw, 'data/processed/counts.bulkized.raw.mtx')
writeMM(counts.bulkized.norm, 'data/processed/counts.bulkized.norm.mtx')

# Make Long ---------------
counts.bulk <- readMM('data/raw/counts.bulk.mtx')

counts.bulk <- counts.bulk %>% as.matrix %>% as.data.frame
counts.bulkized.raw <- counts.bulkized.raw %>% as.matrix %>% as.data.frame
counts.bulkized.norm <- counts.bulkized.norm %>% as.matrix %>% as.data.frame

# Sample ordering is preserved, but they each use different naming schemes. Normalize them...
colnames(counts.bulkized.raw) <- unique(metadata$sample)
rownames(counts.bulkized.raw) <- read.csv('data/processed/genes.bulkized.csv', header = F)$V1
colnames(counts.bulkized.norm) <- colnames(counts.bulkized.raw)
rownames(counts.bulkized.norm) <- rownames(counts.bulkized.raw)

colnames(counts.bulk) <- colnames(counts.bulkized.raw)
rownames(counts.bulk) <- read.csv('data/raw/genes.bulk.csv', header = F)$V1

# Melt the data into long format...
counts.bulk <- counts.bulk %>% rownames_to_column('gene') %>% melt(value.name = 'bulk', id.vars = 'gene', variable.name = 'sample')
counts.bulkized.raw <- counts.bulkized.raw %>% rownames_to_column('gene') %>% melt(value.name = 'bulkized.raw', id.vars = 'gene', variable.name = 'sample')
counts.bulkized.norm <- counts.bulkized.norm %>% rownames_to_column('gene') %>% melt(value.name = 'bulkized.norm', id.vars = 'gene', variable.name = 'sample')

# Finally, merge the data together and save.
merge(counts.bulk, counts.bulkized.raw, by = c('gene', 'sample')) %>%
  merge(counts.bulkized.norm, by = c('gene', 'sample')) %>% saveRDS('data/processed/counts.merged.long.rds')
rm(counts.bulk, counts.bulkized.raw, counts.bulkized.norm, metadata)

# Helper Functions ---------------

#' loadCounts
#' Load the count matrices into the global environment.
#'
#' @param asLog Logical; whether or not to return counts as log counts. Default is F.
#' @param filterGenes Logical; whether or not to only include genes in common between samples. Default is T.
#' @param bulkized Character; the filename of the bulkized counts to read. Default is 'counts.bulkized.norm.mtx'.
#'
#' @return NULL, but adds expr.bulk and expr.bulkized to the global environment.
loadCounts <- function(asLog = F, filterGenes = T, bulkized = 'counts.bulkized.norm.mtx') {
  expr.bulkized <<- readMM(paste0('data/processed/', bulkized)) %>% as('dgCMatrix')
  expr.bulk <<- readMM('data/raw/counts.bulk.mtx') %>% as('dgCMatrix')
  meta <- read.csv('data/raw/meta.singlecell.csv')
  
  colnames(expr.bulkized) <<- unique(meta$sample)
  colnames(expr.bulk) <<- unique(meta$sample)
  rownames(expr.bulkized) <<- as.character(read.csv('data/processed/genes.bulkized.csv', header = F)$V1)
  rownames(expr.bulk) <<- as.character(read.csv('data/raw/genes.bulk.csv', header = F)$V1)
  
  if(filterGenes) {
    genes <- intersect(rownames(expr.bulkized), rownames(expr.bulk))
    expr.bulkized <<- expr.bulkized[genes, ]
    expr.bulk <<- expr.bulk[genes, ]
  }
  
  if(asLog) {
    expr.bulkized <<- log2(1 + expr.bulkized)
    expr.bulk <<- log2(1 + expr.bulk)
  }
}
