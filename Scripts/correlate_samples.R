library(ggplot2)
library(dplyr)
library(limma)
library(reshape2)

# Load our data.
counts <- readRDS('data/processed/counts.merged.long.rds')

# Log ---------------
# Plot the log-count correlations and score correlations
counts %>% ggplot(aes(log2(1 + bulk), log2(1 + bulkized.norm))) + geom_point(alpha = 0.05) +
  xlab('Log Bulk Counts') + ylab('Log Bulkized Counts') + ggtitle('Log Count Correlation') + facet_wrap(~sample, scales = 'free')
ggsave('results/figures/correlations.log.norm.jpg', width = 8.5, height = 8.5, units = 'in')

correlations.log.norm <- lapply(unique(counts$sample), function(sample)
  cor(log2(1 + counts[counts$sample == sample, c('bulk', 'bulkized.norm')]))[1,2]) %>% unlist
correlations.log.norm <- data.frame(sample = unique(counts$sample), correlation = correlations.log.norm)

correlations.log.raw <- lapply(unique(counts$sample), function(sample)
  cor(log2(1 + counts[counts$sample == sample, c('bulk', 'bulkized.raw')]))[1,2]) %>% unlist
correlations.log.raw <- data.frame(sample = unique(counts$sample), correlation = correlations.log.raw)

# Voom + Quantile Norm ---------------
# Voom and quantile normalize.
counts$bulkized.raw <- counts[, c(1, 2, 4)] %>% # Drop the bulk expression
  dcast(gene~sample, sum, value.var = 'bulkized.raw') %>% # Make the data wide
  .[, -1] %>% # Drop gene IDs
  voom(normalize.method = 'quantile') %>% .[[1]] %>% # Voom/quantile normalize
  melt(value.name = 'bulkized.raw') %>% .[, 'bulkized.raw'] # Melt again and extract the bulkized values.

counts$bulkized.norm <- counts[, c(1, 2, 5)] %>% # Drop the bulk expression
  dcast(gene~sample, sum, value.var = 'bulkized.norm') %>% # Make the data wide
  .[, -1] %>% # Drop gene IDs
  voom(normalize.method = 'quantile') %>% .[[1]] %>% # Voom/quantile normalize
  melt(value.name = 'bulkized.norm') %>% .[, 'bulkized.norm'] # Melt again and extract the bulkized values.

counts$bulk <- counts[, 1:3] %>% # Drop the bulkized expression
  dcast(gene~sample, sum, value.var = 'bulk') %>% # Make the data wide
  .[, -1] %>% # Drop gene IDs
  voom(normalize.method = 'quantile') %>% .[[1]] %>% # Voom/quantile normalize
  melt(value.name = 'bulk') %>% .[, 'bulk'] # Melt again and extract the bulkized values.

# Plot the voomed-count correlations and score correlations
counts %>% ggplot(aes(bulk, bulkized.norm)) + geom_point(alpha = 0.05) +
  xlab('Normalized Bulk Counts') + ylab('Normalized Bulkized Counts') + ggtitle('Normalized Count Correlation') + facet_wrap(~sample, scales = 'free')
ggsave('results/figures/correlations.log_normal.norm.jpg', width = 8.5, height = 8.5, units = 'in')

correlations.voom.norm <- lapply(unique(counts$sample), function(sample)
  cor(counts[counts$sample == sample, c('bulk', 'bulkized.norm')])[1,2]) %>% unlist()
correlations.voom.norm <- data.frame(sample = unique(counts$sample), correlation = correlations.voom.norm)

correlations.voom.raw <- lapply(unique(counts$sample), function(sample)
  cor(counts[counts$sample == sample, c('bulk', 'bulkized.raw')])[1,2]) %>% unlist()
correlations.voom.raw <- data.frame(sample = unique(counts$sample), correlation = correlations.voom.raw)

# Save ---------------
write.table(correlations.log.norm, 'results/correlation.log.norm.csv', row.names = F, quote = F, sep = ', ')
write.table(correlations.log.raw, 'results/correlation.log.raw.csv', row.names = F, quote = F, sep = ', ')
write.table(correlations.voom.norm, 'results/correlation.log_normal.norm.csv', row.names = F, quote = F, sep = ', ')
write.table(correlations.voom.raw, 'results/correlation.log_normal.raw.csv', row.names = F, quote = F, sep = ', ')

rm(correlations.log.raw, correlations.log.norm, correlations.voom.raw, correlations.voom.norm, counts)