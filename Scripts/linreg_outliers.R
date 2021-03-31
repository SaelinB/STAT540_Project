library(dplyr)
library(reshape2)
library(tibble)
library(Matrix)
library(tidyverse)
library(pracma)
library(data.table)
library(NMF)
library(parallel)
options('mc.cores' = 4) # Use 4 cores

# Load our data
loadCounts(bulkized="counts.bulkized.norm.rds")

# Fit a regression line along bulk~bulkized
genes.residuals <- do.call(cbind, lapply(colnames(expr.bulk), function(sample) {
  expr.bulk.log <- log2(1 + expr.bulk[, sample])
  expr.bulkized.log <- log2(1 + expr.bulkized[, sample])
  lm.summary <- lm(expr.bulk.log~expr.bulkized.log) %>% summary
  
  lm.summary$residuals
}))
colnames(genes.residuals) <- colnames(expr.bulkized)
rownames(genes.residuals) <- rownames(expr.bulkized)

# Plot the histogram of residuals faceted by sample
genes.residuals.long <- reshape2::melt(genes.residuals, measure.vars=colnames(expr.bulkized))
colnames(genes.residuals.long) <- c("gene", "sample", "res")
genes.residuals.long %>% ggplot(aes(res)) + geom_histogram(bins = 60) + facet_wrap(~sample, scales = 'free_y') + xlim(-2, 2)

# Compile a list of outliers (<1st quartile and >3rd quartile) for every sample and among all samples
genes.sample_outliers.first_quartile <- list()
genes.sample_outliers.third_quartile <- list()
genes.sample_outliers.both_quartiles <- list()
genes.all_outliers.first_quartile <- list()
genes.all_outliers.third_quartile <- list()
genes.all_outliers.both_quartiles <- list()
invisible(lapply(colnames(expr.bulkized), function(sample) {
  ordered_residuals <- data.frame(genes.residuals[, sample])
  ordered_residuals <- cbind(rownames(ordered_residuals), ordered_residuals)
  colnames(ordered_residuals) <- c("gene", "res")
  rownames(ordered_residuals) <- NULL
  ordered_residuals <- data.frame(ordered_residuals[order(ordered_residuals$res),])
  
  threshold_left <- quantile(ordered_residuals$res, 0.25)
  threshold_right <- quantile(ordered_residuals$res, 0.75)
  
  first_quartile_residuals <- ordered_residuals[ordered_residuals$res < threshold_left,]
  third_quartile_residuals <- ordered_residuals[ordered_residuals$res > threshold_right,]
  both_quartiles_residuals <- ordered_residuals[ordered_residuals$res < threshold_left | ordered_residuals$res > threshold_right,]
  
  genes.sample_outliers.first_quartile[[sample]] <<- first_quartile_residuals$gene %>% as.character()
  genes.sample_outliers.third_quartile[[sample]] <<- third_quartile_residuals$gene %>% as.character()
  genes.sample_outliers.both_quartiles[[sample]] <<- both_quartiles_residuals$gene %>% as.character()
  
  genes.all_outliers.first_quartile <<- union(genes.all_outliers.first_quartile, genes.sample_outliers.first_quartile[[sample]])
  genes.all_outliers.third_quartile <<- union(genes.all_outliers.third_quartile, genes.sample_outliers.third_quartile[[sample]])
  genes.all_outliers.both_quartiles <<- union(genes.all_outliers.both_quartiles, genes.sample_outliers.both_quartiles[[sample]])
}))

# Count how many samples outliers are found in
genes.common_outliers.sample_appearances.first_quartile <- do.call(rbind, lapply(names(genes.sample_outliers.first_quartile), function(sample) data.frame(sample = sample, gene = genes.sample_outliers.first_quartile[[sample]]))) %>%
  as.data.table %>% .[, .N, gene]
genes.common_outliers.sample_appearances.third_quartile <- do.call(rbind, lapply(names(genes.sample_outliers.third_quartile), function(sample) data.frame(sample = sample, gene = genes.sample_outliers.third_quartile[[sample]]))) %>%
  as.data.table %>% .[, .N, gene]
genes.common_outliers.sample_appearances.both_quartiles <- do.call(rbind, lapply(names(genes.sample_outliers.both_quartiles), function(sample) data.frame(sample = sample, gene = genes.sample_outliers.both_quartiles[[sample]]))) %>%
  as.data.table %>% .[, .N, gene]

# Make a list of outliers found in 2/3 of samples
genes.common_outliers.first_quartile <- genes.common_outliers.sample_appearances.first_quartile[N >= length(colnames(expr.bulkized)) * 2 / 3]$gene
genes.common_outliers.third_quartile <- genes.common_outliers.sample_appearances.third_quartile[N >= length(colnames(expr.bulkized)) * 2 / 3]$gene
genes.common_outliers.both_quartiles <- genes.common_outliers.sample_appearances.both_quartiles[N >= length(colnames(expr.bulkized)) * 2 / 3]$gene

print("First quartile common and all outliers:")
length(genes.common_outliers.first_quartile)
length(genes.all_outliers.first_quartile)

print("Third quartile common and all outliers:")
length(genes.common_outliers.third_quartile)
length(genes.all_outliers.third_quartile)

print("Both quartiles common and all outliers:")
length(genes.common_outliers.both_quartiles)
length(genes.all_outliers.both_quartiles)

# Load cell types' marker genes and intersect it with actually expressed genes
genes.expressed <- rownames(expr.bulkized)
marker.genes <- readRDS('data/processed/genes.markers.rds')
marker.genes.common <- list()
invisible(lapply(names(marker.genes), function(cell_type) {
  marker_genes <- unique(marker.genes[[cell_type]]$ensembl_gene_id)
  marker.genes.common[[cell_type]] <<- intersect(marker_genes, genes.expressed)
}))

# Make a dataframe with cell types, corresponding marker genes, samples and their residuals
marker.genes.cell_type <- do.call(rbind, lapply(names(marker.genes.common), function(cell_type) {
  df <- data.frame(`cell type` = cell_type, gene = marker.genes.common[[cell_type]])
  df <- df[order(df$gene),]
  res <- genes.residuals[rownames(genes.residuals) %in% marker.genes.common[[cell_type]],]
  res <- cbind(rownames(res), res)
  res <- res[order(res[, 1]),]
  cbind(df, res[, -1])
}))
marker.genes.cell_type <- reshape2::melt(marker.genes.cell_type, measure.vars=colnames(expr.bulkized))
colnames(marker.genes.cell_type) <- c("cell_type", "gene", "sample", "res")
marker.genes.cell_type[, "res"] <- sapply(marker.genes.cell_type[, "res"], as.numeric)

# Sort genes by residual values
marker.genes.cell_type <- marker.genes.cell_type[order(marker.genes.cell_type$res),]
genes.residuals.long <- genes.residuals.long[order(genes.residuals.long$res),]

# Make genes factors and reorder factor values by residuals for nice plots
marker.genes.cell_type$gene <- factor(marker.genes.cell_type$gene) %>% fct_reorder(marker.genes.cell_type$res)
genes.residuals.long$gene <- factor(genes.residuals.long$gene) %>% fct_reorder(genes.residuals.long$res)

# Plot residuals of all genes for all samples
genes.residuals.long %>% ggplot(aes(x=gene, y=res, fill=sample)) +
  geom_col() +
  ylab("Residual") +
  xlab("Genes") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) #+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('results/figures/residuals.all_genes.png', width = 8.5, height = 8.5, units = 'in')

# Plot residuals of marker genes faceted by cell type
marker.genes.cell_type %>% ggplot(aes(x=gene, y=res, fill=sample)) +
  geom_col() +
  ylab("Residual") +
  xlab("Marker genes") +
  facet_wrap(~cell_type, scales="free_x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) #+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('results/figures/residuals.marker_genes.png', width = 8.5, height = 8.5, units = 'in')

# For each cell type, get the number of marker genes that are outliers (separate matrices for down and upregulated)
bulkized.upregulated <- do.call(rbind, lapply(names(marker.genes.common), function(cell_type) {
  marker_genes <- marker.genes.common[[cell_type]]
  data.frame(`cell type` = cell_type, `marker genes` = length(intersect(marker_genes, genes.common_outliers.first_quartile)), `marker genes total` = length(marker_genes), percent = length(intersect(marker_genes, genes.common_outliers.first_quartile)) * 100 / length(marker_genes))
})) %>% setorder(-percent)
bulkized.downregulated <- do.call(rbind, lapply(names(marker.genes.common), function(cell_type) {
  marker_genes <- marker.genes.common[[cell_type]]
  data.frame(`cell type` = cell_type, `marker genes` = length(intersect(marker_genes, genes.common_outliers.third_quartile)), `marker genes total` = length(marker_genes), percent = length(intersect(marker_genes, genes.common_outliers.third_quartile)) * 100 / length(marker_genes))
})) %>% setorder(-percent)

# Make a list of number of marker genes for each cell type
markers.N <- lapply(names(marker.genes), function(cell_type) {
  length(unique(marker.genes[[cell_type]]$ensembl_gene_id))
})
names(markers.N) <- names(marker.genes)
marker.names <- names(markers.N)
names(marker.names) <- names(markers.N)

# Sort cell type outlier marker genes by cell type name
bulkized.downregulated.sorted <- bulkized.downregulated[, 1:2]
bulkized.downregulated.sorted <- bulkized.downregulated.sorted[order(bulkized.downregulated.sorted$cell.type),]
bulkized.upregulated.sorted <- bulkized.upregulated[, 1:2]
bulkized.upregulated.sorted <- bulkized.upregulated.sorted[order(bulkized.upregulated.sorted$cell.type),]

singlecell.metadata <- fread('data/raw/meta.singlecell.csv')

marker.map <- list(Astrocyte = c('AST-FB', 'AST-PP'),
                   Endothelial = 'Endothelial',
                   GabaPV = 'IN-PV',
                   GabaRelnCalb = '',
                   GabaVIPReln = 'IN-VIP',
                   Layer.2.3.Pyra = 'L2/3',
                   Layer.4.Pyra = 'L4',
                   Layer.6a.Pyra = '',
                   Layer.6b.Pyra = '',
                   Microglia = 'Microglia',
                   Microglia_activation = 'Microglia',
                   Microglia_deactivation = 'Microglia',
                   Oligo = 'Oligodendrocytes',
                   OligoPrecursors = 'OPC',
                   Pyramidal_Glt_25d2 = '', # Layer 5
                   Pyramidal_S100a10 = '', # Layer 5
                   PyramidalCorticoThalam = '', # Layer 6
                   Pyramidal = c('L5/6', 'L5/6-CC'))

cell_type.names <- list(Astrocyte = 'Astrocyte',
                        Endothelial = 'Endothelial',
                        GabaPV = 'Fast Spiking Basket',
                        GabaRelnCalb = 'Martinotti',
                        GabaVIPReln = 'GABAergic VIP Reelin',
                        Layer.2.3.Pyra = 'L2/3',
                        Layer.4.Pyra = 'L4',
                        Layer.6a.Pyra = 'L6a',
                        Layer.6b.Pyra = 'L6b',
                        Microglia = 'Microglia',
                        Microglia_activation = 'Microglia (Activated)',
                        Microglia_deactivation = 'Microglia (Deactivated)',
                        Oligo = 'Oligodendrocyte',
                        OligoPrecursors = 'OPC',
                        Pyramidal_Glt_25d2 = 'Pyramidal (GLT25D2)',
                        Pyramidal_S100a10 = 'Pyramidal (P11)',
                        PyramidalCorticoThalam = 'Pyramidal (cortico-thalamic)',
                        Pyramidal = 'Pyramidal')

# Combine the lists for up and downregulated markers for each cell type
up.down.info <- cbind(bulkized.upregulated.sorted, bulkized.downregulated.sorted[, 2])
up.down.info <- t(up.down.info)
colnames(up.down.info) <- up.down.info[1, ]
up.down.info <- data.frame(up.down.info[2:3, ])

cell_count <- lapply(colnames(up.down.info), function(cell_type) {
 singlecell.metadata[cluster %in% marker.map[[cell_type]], .N]
}) %>% data.frame()
colnames(cell_count) <- colnames(up.down.info)
up.down.info <- rbind(cell_count, up.down.info)

up.down.info <- sapply(up.down.info, strtoi)

rownames(up.down.info) <- c("Cells", "Up", "Down")
colnames(up.down.info) <- lapply(colnames(up.down.info), function(cell_type) {
  cell_type.names[[cell_type]]  
})

# Correlate num of outlier marker genes with the cell count of that cell type
cor.test(up.down.info[1,], up.down.info[3,])

# TODO Run the first line and make the following edits to make the colnames rotate 45 degrees.
# trace(NMF:::draw_colnames, edit = T)
# rot <- 45
# vjust <- 0.5
# hjust <- 1
# y <- unit(1, "npc") - unit(10, "bigpts")
aheatmap(t(t(up.down.info) / unlist(markers.N)) * 100,
         breaks = seq(0, max(t(t(up.down.info[-1, ]) / unlist(markers.N)) * 100), length.out = 10),
         border_color = 'black', txt = up.down.info,
         color = 'Blues', Rowv = F,
         Colv = order(-t(t(up.down.info) / unlist(markers.N))[3,]),
         labCol = unlist(cell_type.names), main = 'Cell Type Markers',
         cellwidth = 35, cellheight = 35, width = 11, filename = 'results/figures/linreg_underexpressed_celltypes.png')

# Get marker genes for each cell type that are outliers in both directions
bulkized.upanddownregulated <- do.call(rbind, lapply(names(marker.genes.common), function(cell_type) {
  marker_genes <- marker.genes.common[[cell_type]]
  data.frame(`cell type` = cell_type, `marker genes` = length(intersect(marker_genes, genes.common_outliers.both_quartiles)), `marker genes total` = length(marker_genes), percent = length(intersect(marker_genes, genes.common_outliers.both_quartiles)) * 100 / length(marker_genes))
})) %>% setorder(-percent)

# Mann Whitney test to compare the distribution of residuals for all genes and for each set of marker genes
mann_whitney <- do.call(rbind, lapply(names(marker.names), function(x) {
  data.frame(cell_type = x, p.value = wilcox.test(genes.residuals.long$res, marker.genes.cell_type[marker.genes.cell_type$cell_type == x,]$res)$p.value, alternative = 'greater')
}))
mann_whitney <- cbind(mann_whitney, data.frame(adj.p.value = p.adjust(mann_whitney$p.value, 'BH', length(names(marker.names)))))
mann_whitney

print(paste0("% of outlier downregulated genes that are markers: ", (sum(bulkized.downregulated.sorted$marker.genes) * 100 / length(genes.common_outliers.third_quartile)), "%"))
print(paste0("% of marker genes that are not downregulated outliers: ", ((Reduce("+", markers.N) - sum(bulkized.downregulated.sorted$marker.genes)) * 100 / Reduce("+", markers.N)), "%"))
print(paste0("% of marker genes that are not upregulated outliers: ", ((Reduce("+", markers.N) - sum(bulkized.upregulated.sorted$marker.genes)) * 100 / Reduce("+", markers.N)), "%"))
print(paste0("% of marker genes that are not at all outliers: ",
  ((Reduce("+", markers.N) - sum(bulkized.downregulated.sorted$marker.genes)  - sum(bulkized.upregulated.sorted$marker.genes)) * 100 / Reduce("+", markers.N)), "%"))
