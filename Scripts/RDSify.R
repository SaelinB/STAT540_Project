library(dplyr)
library(Matrix)

# Simply load counts...
loadCounts()

# And save them with a different format!
saveRDS(expr.bulk %>% as.matrix, 'data/raw/counts.bulk.rds')
saveRDS(expr.bulkized %>% as.matrix, 'data/processed/counts.bulkized.norm.rds')
loadCounts(bulkized = 'counts.bulkized.raw.mtx')
saveRDS(expr.bulkized %>% as.matrix, 'data/processed/counts.bulkized.raw.rds')

# Also re-source the loadCounts function to work with the RDS objects instead
loadCounts <- function(asLog = F, bulkized = 'counts.bulkized.norm.rds') {
  expr.bulkized <<- readRDS(paste0('data/processed/', bulkized))
  expr.bulk <<- readRDS('data/raw/counts.bulk.rds')
  meta <- read.csv('data/raw/meta.singlecell.csv')
  
  if(asLog) {
    expr.bulkized <<- log2(1 + expr.bulkized)
    expr.bulk <<- log2(1 + expr.bulk)
  }
}
