library(dplyr)
library(reshape2)
library(tibble)
library(Matrix)

# Load our data ---------------
loadCounts()

# Within-sample corrections ---------------
expr.bulkized.corrected.log <- do.call(cbind, lapply(colnames(expr.bulk), function(sample) {
  expr.bulk.log <- log2(1 + expr.bulk[, sample])
  expr.bulkized.log <- log2(1 + expr.bulkized[, sample])
  lm.summary <- lm(expr.bulk.log~expr.bulkized.log) %>% summary
  
  expr.bulkized.log * lm.summary$coefficients[2, 1] +
    lm.summary$coefficients[1, 1]
}))

colnames(expr.bulkized.corrected.log) <- colnames(expr.bulkized)
rownames(expr.bulkized.corrected.log) <- rownames(expr.bulkized)

# Between-sample corrections ---------------
expr.bulkized.corrected.log.means <- Matrix(as.matrix(expr.bulkized.corrected.log)) %>% Matrix::rowMeans()
expr.bulk.log.means <- Matrix(as.matrix(log2(1 + expr.bulk))) %>% Matrix::rowMeans()

lm.summary <- lm(expr.bulk.log.means~expr.bulkized.corrected.log.means) %>% summary
expr.bulkized.corrected.log.normalized <- expr.bulkized.corrected.log * lm.summary$coefficients[2, 1] +
  lm.summary$coefficients[1, 1]

expr.bulkized.corrected <- pmax(2 ^ expr.bulkized.corrected.log - 1, 0) %>% Matrix()
expr.bulkized.corrected.between <- pmax(2 ^ expr.bulkized.corrected.log.normalized - 1, 0) %>% Matrix()

# Save ---------------
saveRDS(expr.bulkized.corrected, 'data/processed/counts.bulkized.corrected.rds')
saveRDS(expr.bulkized.corrected.between, 'data/processed/counts.bulkized.corrected.between.rds')
