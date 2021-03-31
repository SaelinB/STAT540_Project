library(dplyr)
library(limma)
library(Matrix)
library(data.table)

# Load counts ---------------
loadCounts(asLog = T)

# Scale so the means are equal per sample. This reduces undesirable data type variance
expr.bulkized.scaled <- t(t(expr.bulkized) * colMeans(expr.bulk) / colMeans(expr.bulkized))

# Bind them together for DE analysis
expr.merged <- cbind(expr.bulk, expr.bulkized)
rm(expr.bulk, expr.bulkized)

# DE Analysis ---------------
# Setup the covariates
style <- factor(c(rep('Bulk', 41), rep('Bulkized', 41)), levels = c('Bulk', 'Bulkized'))
sample <- colnames(expr.merged) %>% as.factor
region <- read.csv('data/raw/meta.bulk.csv')$region %>% as.factor

# Do DE analysis between styles (bulk vs. bulkized) using individual as a covariate
# Negative logFC means "higher in bulk"
DEGs <- lmFit(expr.merged, model.matrix(~style + sample + region)) %>% eBayes %>%
  topTable(coef = 'styleBulkized', number = Inf, sort.by = 'P')

# Plot the p-value distributions
DEGs %>% ggplot(aes(P.Value)) +
  geom_histogram(bins = 60, fill = 'white', color = 'black') +
  xlab('P-Value') + ylab('Count') + ggtitle('P-Value Distribution') +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1000))
ggsave('results/figures/p.distribution.pdf', width = 8.5, height = 8.5, units = 'in')

# Save ---------------
saveRDS(DEGs, 'data/processed/DEGs.rds')

rm(DEGs, expr.merged, style, sample, region)
