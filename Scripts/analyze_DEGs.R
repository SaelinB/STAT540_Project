library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tibble)
library(ermineR)
library(NMF)

# Load data ---------------
metadata <- fread('data/raw/meta.singlecell.csv')
DEGs <- readRDS('data/processed/DEGs.rds')
DEGs$primerid <- rownames(DEGs)
markers <- readRDS('data/processed/genes.markers.rds') %>% lapply(as.data.table)
DEGs <- DEGs %>% as.data.table

# Make a name map so we can pretty print
marker.names <- list(Astrocyte = 'Astrocyte',
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

# Subset for genes in our data
genes.dataset <- readRDS('data/processed/counts.merged.long.rds')$gene %>% unique
markers.all <- lapply(markers, function(marker) marker[, ensembl_gene_id]) %>% unlist %>% unique
markers <- lapply(markers, function(marker) marker[ensembl_gene_id %in% genes.dataset])
markers.N <- lapply(markers, function(marker) length(unique(marker$ensembl_gene_id)))

# Only maintain significant DEs with a logFC of at least a quarter of the maximal
DEGs.filtered <- DEGs %>% .[adj.P.Val < 0.05 & abs(logFC) > (0.25 * max(abs(logFC)))]

# Print information on the number that are down-regulated and up-regulated
DEGs.filtered[primerid %in% markers.all, .N, logFC < 0]

# Make a name map so we can associate cell types from NeuroExpresso with cell types annotated in our data
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

# Compute markers that are up and down regulated, and associate it with their cell type.
up.down.info <- do.call(rbind, lapply(names(markers), function(marker) {
  up <- DEGs.filtered[primerid %in% markers[[marker]]$ensembl_gene_id & logFC > 0, primerid] %>% length
  down <- DEGs.filtered[primerid %in% markers[[marker]]$ensembl_gene_id & logFC < 0, primerid] %>% length
  
  if(length(up) == 0) up <- 0
  if(length(down) == 0) down <- 0
  
  data.frame(Cell.type = marker, Cells = metadata[cluster %in% marker.map[[marker]], .N], Up = up, Down = down)
})) %>% column_to_rownames('Cell.type') %>% as.matrix %>% t

# Assess correlation between number of cells per cell type and number of down-regulated
cor.test(up.down.info[1, ], up.down.info[3, ])

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
         labCol = unlist(marker.names), main = 'Cell Type Markers',
         cellwidth = 35, cellheight = 35, width = 11, filename = 'results/figures/DEGs.pdf')

# Tabulate counts per cell type ---------------
# Marker counts
lapply(markers, function(set) {
  DEGs.filtered[logFC < 0 & primerid %in% set[, unique(ensembl_gene_id)], set[ensembl_gene_id %in% primerid]]
}) -> marker.populations

# Background counts
lapply(markers, function(set) {
  DEGs.significant[primerid %in% set[, unique(ensembl_gene_id)], set[ensembl_gene_id %in% primerid]]
}) -> background.populations

marker.counts <- lapply(marker.populations, function(set) {
  list(N = length(unique(set[, ensembl_gene_id])),
       GO = set[, .N, go_id],
       chromosome = unique(set[, chromosome_name, ensembl_gene_id])[, .N, chromosome_name],
       length = unique(set[, transcript_length, ensembl_gene_id])[, transcript_length],
       GC = unique(set[, percentage_gene_gc_content, ensembl_gene_id])[, percentage_gene_gc_content])
})

background.counts <- lapply(background.populations, function(set) {
  list(N = length(unique(set[, ensembl_gene_id])),
       GO = set[, .N, go_id],
       chromosome = unique(set[, chromosome_name, ensembl_gene_id])[, .N, chromosome_name],
       length = unique(set[, transcript_length, ensembl_gene_id])[, transcript_length],
       GC = unique(set[, percentage_gene_gc_content, ensembl_gene_id])[, percentage_gene_gc_content])
})

saveRDS(marker.populations, 'data/processed/DEGs.markers.rds')
saveRDS(marker.counts, 'data/processed/DEGs.markers.counts.rds')

theme_set(theme_cowplot())

# Plot the number of DEGs by cell type, colored by cell type
do.call(rbind, lapply(marker.counts, function(x) c(N = x$N))) %>% cbind(marker.N = unlist(markers.N)) %>%
  cbind(do.call(rbind, lapply(background.counts, function(x) c(bg.N = x$N)))) %>%
  as.data.frame %>% rownames_to_column('CellType') %>% as.data.table %>%
  .[, Percent := 100 * N / marker.N / bg.N] %>%
  ggplot(aes(CellType, Percent, fill = CellType)) + geom_bar(stat = 'identity') + coord_flip() +
  xlab('Cell Type') + ggtitle('DEGs by Cell Type') + theme(legend.position = 'none') +
  scale_y_continuous(expand = c(0, 0))
ggsave('data/processed/figures/DEGs.CellType.pdf', width = 8.5, height = 11, units = 'in')

# Plot the number of DEGs by GC content, colored by cell type
do.call(rbind, lapply(names(marker.counts), function(x) {
  if(length(marker.counts[[x]]$length) != 0)
    data.frame(CellType = x, gc = marker.counts[[x]]$GC)
})) %>% as.data.frame %>% ggplot(aes(x = gc, fill = CellType)) + geom_density(alpha = 0.1) +
  ggtitle('DEGs by GC Content') + theme(legend.position = c(0.7, 0.7)) + ylab('N') + xlab('GC (%)') +
  scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
ggsave('data/processed/figures/DEGs.GC.pdf', width = 11, height = 8.5, units = 'in')

# Plot the number of DEGs by the length of the gene, colored by cell type
do.call(rbind, lapply(names(marker.counts), function(x) {
  if(length(marker.counts[[x]]$length) != 0)
    data.frame(CellType = x, length = marker.counts[[x]]$length)
})) %>% as.data.frame %>% ggplot(aes(length, fill = CellType)) + geom_histogram() +
  ggtitle('DEGs by Length') + theme(legend.position = c(0.7, 0.7)) + ylab('N') + xlab('Length (bp)') +
  scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0), limits = c(0, 6000))
ggsave('data/processed/figures/DEGs.Length.pdf', width = 11, height = 8.5, units = 'in')

# Plot the number of DEGs by the chromosome they're on, colored by cell type
do.call(rbind, lapply(names(marker.counts), function(x) {
  if(length(marker.counts[[x]]$length) != 0)
    data.frame(CellType = x, chromosome = marker.counts[[x]]$chromosome$chromosome_name, N = marker.counts[[x]]$chromosome$N)
})) %>% as.data.frame() %>% ggplot(aes(chromosome, N, fill = CellType)) + geom_bar(stat = 'identity') +
  theme(legend.position = 'none') + coord_flip() + xlab('Chromosome') + ggtitle('DEGs by Chromosome') +
  scale_y_continuous(expand = c(0, 0))
ggsave('data/processed/figures/DEGs.Chromosome.pdf', width = 11, height = 8.5, units = 'in')

# Gene Enrichment Analysis ---------------
genes <- readRDS('data/processed/genes.detected.rds') %>% as.data.table
goi <- genes[ensembl_gene_id %in% DEGs.filtered$primerid, .(ensembl_gene_id, hgnc_symbol)] %>%
  merge(DEGs.filtered, by.x = 'ensembl_gene_id', by.y = 'primerid') %>% unique %>%
  .[hgnc_symbol != '' & hgnc_symbol != 'PINX1'] %>% as.data.frame %>%
  column_to_rownames('hgnc_symbol')
goi$logFC <- abs(goi$logFC)

oraOut <- ora(annotation = 'Generic_human',
              scores = goi,
              scoreColumn = 'logFC',
              bigIsBetter = T)
