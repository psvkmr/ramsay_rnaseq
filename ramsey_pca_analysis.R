library(PCAtools)
library(DESeq2)
library(tidyverse)
library(vsn)

#https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html

setwd('ramsey/')

# load count matrices
load('ramsey_s303_norm_counts_matrix.RData')
load('ramsey_s303_non_zero_counts_matrix.RData')

# load metadata
md <- read.csv('s303_metadata.csv')

# add matrix sample names to metadata variable as different
#md$matrix_id <- stringr::str_sort(colnames(zcm), numeric = T)
# wrong length because NAT 4 is in metadata but not matrix
md <- md[md$sample_id != 'NAT 4', ]
# works after removing
md$matrix_id <- stringr::str_sort(colnames(zcm), numeric = T)
# check to make sure
#View(md[, c('sample_id', 'matrix_id')])

# re-sort whole md to match matrix variable order
md <- arrange(md, factor(matrix_id, levels = colnames(zcm)))


# convert matrix ids to row names
md <- column_to_rownames(md, 'matrix_id')

# edit md 'reactive lymph nodes' to shorten
md$sample_group <- gsub('Lymph Nodes\n$', 'LN', md$sample_group)
md$sample_group <- gsub('\n$', '', md$sample_group)

# in this analysis, only the 'NAT' samples are going to be used, 
# so both the md and cm will be subsetted for these
smpls <- colnames(zcm)[grep('NAT', colnames(zcm))]
nat.cm <- zcm[, smpls]
nat.md <- md[smpls, ]


# sample groups
smpl.grps <- list(cll_ln = rownames(nat.md[nat.md$sample_group == 'CLL LN', ]), 
                  cll_pb = rownames(nat.md[nat.md$sample_group == 'CLL PB', ]), 
                  ctrl = rownames(nat.md[nat.md$sample_group == 'Control', ]), 
                  reac = rownames(nat.md[nat.md$sample_group == 'Reactive LN', ]))

# cms subset by sample group comparisons
zcms <- list(ctrl_v_cll_ln = zcm[, c(smpl.grps$ctrl, smpl.grps$cll_ln)], 
             ctrl_v_cll_pb = zcm[, c(smpl.grps$ctrl, smpl.grps$cll_pb)], 
             ctrl_v_reac = zcm[, c(smpl.grps$ctrl, smpl.grps$reac)])
#lapply(zcms, dim)

# cms filtered for low counts specifically in comparison subset
zcms <- lapply(zcms, function(mat) mat[rowSums(mat) > (((dim(mat)[2])/2)-1), ])

# normalised subset cms
ncms <- list(ctrl_v_cll_ln = ncm[, c(smpl.grps$ctrl, smpl.grps$cll_ln)], 
             ctrl_v_cll_pb = ncm[, c(smpl.grps$ctrl, smpl.grps$cll_pb)], 
             ctrl_v_reac = ncm[, c(smpl.grps$ctrl, smpl.grps$reac)])
#lapply(ncms, dim)


# outlier check
outlierCheck <- function(mat){
  g <- list()
  s <- list()
  for (i in 1:(dim(mat)[1])){
    mean.i <- mean(mat[i, ])
    sd.i <- sd(mat[i, ])
    for (j in mat[i, ]){
      if (j > (mean.i + (3*sd.i)) & sd.i > 1){
        g <- append(g, row.names(mat)[i])
        s <- append(s, colnames(mat)[which(mat[i, ] == j)])
      } else if (j < (mean.i - (3*sd.i)) & sd.i > 1){
        g <- append(g, row.names(mat)[i])
        s <- append(s, colnames(mat)[which(mat[i, ] == j)])
      }
    }
  }
  o <- data.frame('sample' = unlist(s), 'gene' = unlist(g))
  return(o)
}

outliers <- lapply(ncms, outlierCheck)
#table(outliers$ctrl_v_cll_pb$sample )
# suggests potential outlier sample NAT.16_S16

plotOutliers <- function(outlier){
  df <- as.data.frame(table(outlier$sample))
  if (nrow(df) != 0){ 
    ggplot(df, aes(Var1, Freq)) + 
      geom_col() +
      labs(x = 'Sample', y = 'Count') + 
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  }
}

outliers.plots <- lapply(outliers, plotOutliers)

# subset metadata by comparisons
mds <- list(ctrl_v_cll_ln = md[c(smpl.grps$ctrl, smpl.grps$cll_ln), ], 
            ctrl_v_cll_pb = md[c(smpl.grps$ctrl, smpl.grps$cll_pb), ], 
            ctrl_v_reac = md[c(smpl.grps$ctrl, smpl.grps$reac), ])

#lapply(mds, dim)

# stabilise variance across data to make it homoskedastic
svcms <- lapply(zcms, varianceStabilizingTransformation)
# similar to scale(log2(ncms$ctrl_v_reac + 1), center = T, scale = T)

# check variance transformation adequate
sdPlots <- lapply(svcms, meanSdPlot)

# plot ranges of log counts per sample for matrix subsets
boxplot(log2(svcms$ctrl_v_cll_ln))
boxplot(log2(svcms$ctrl_v_cll_pb))
boxplot(log2(svcms$ctrl_v_reac))

# compare total vst counts between two samples
plot(x = svcms$ctrl_v_reac[, 'NAT.5_S5'], y = svcms$ctrl_v_reac[, 'NAT.8_S8'])

# run pca, remove bottom 10% of variables based on variance
pcas <- list()
for (name in names(svcms)){
  p <- pca(svcms[[name]], metadata = mds[[name]], removeVar = 0.1)
  p <- list(p)
  names(p) <- name
  pcas <- append(pcas, p)
}

# plot the components eigenvalues as % variance explanation
screes <- lapply(pcas, function(p) screeplot(p, getComponents(p, 1:8)))

# biplot of all samples by PC1 vs PC2, with clouds
bi.plt.smpls <- 
  lapply(pcas, function(p) biplot(p, showLoadings = F, colby = 'sample_group', encircle = T,
                                  widthConnectors = 0.1, gridlines.minor = F,  gridlines.major = F,
                                  legendPosition = 'right', legendLabSize = 8, legendTitleSize = 8, 
                                  legendIconSize = 4, axisLabSize = 8))


# add gene labels with biggest effects
bi.plt.load <- 
  lapply(pcas, function(p) biplot(p, showLoadings = T, sizeLoadingsNames = 3, lab = NULL, colby = 'sample_group', 
                                  widthConnectors = 0.1, gridlines.minor = F, gridlines.major = F, 
                                  legendPosition = 'right', legendLabSize = 8, legendTitleSize = 8, legendIconSize = 4, 
                                  axisLabSize = 10))

# comparison scatters for PC1-4 on grid
pairs.plt <- 
  lapply(pcas, function(p) pairsplot(p, components = getComponents(p, 1:4), triangle = F, 
                                     colby = 'sample_group', axisLabSize = 8))


# plot individual contribution of genes with largest loadings
loadings.plt <- 
  lapply(pcas, function(p) plotloadings(p, components = getComponents(p, 1:4), rangeRetain = 0.001, 
                                        shapeSizeRange = c(2, 10), absolute = T, legendPosition = 'none', 
                                        #legendLabSize = 8, legendIconSize = 4, 
                                        axisLabSize = 8, gridlines.major = F, gridlines.minor = F))

# plot loadings of each variable for each component
# unnamed columns
eigencor.plt <- 
  lapply(pcas, function(p) eigencorplot(p, getComponents(p, 1:5),
                                        metavars = c("sample_group", "rin_2", "rna_conc_qubit", "rna_dil_conc_qubit", "total_rna_qubit",
                                                     "rna_conc_ndrop", "total_rna_ndrop", "a260_280_ndrop", "a260_230_ndrop"), 
                                        col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'), 
                                        colCorval = 'white', posColKey = 'top'))

# get pearson r2 values and significane of correlation of PCs with metadata variables
# unnamed columns
eigencorr2.plt <-
  lapply(pcas, function(p) eigencorplot(p, getComponents(p, 1:5),
                                        metavars = c('sample_group', "rin_2", "rna_conc_qubit", "rna_dil_conc_qubit", "total_rna_qubit",
                                                     "rna_conc_ndrop", "total_rna_ndrop", "a260_280_ndrop", "a260_230_ndrop"), 
                                        colCorval = 'white', posColKey = 'top', 
                                        plotRsquared = T, corFUN = 'pearson', corUSE = 'pairwise.complete.obs', 
                                        corMultipleTestCorrection = 'BH'))

  
# sample contribution to each PC
rotateds <- lapply(pcas, `[[`, 'rotated')
# also suggests NAT7 outlier

# each gene contribution to each PC
loads <- lapply(pcas, `[[`, 'loadings')
# as.data.frame(loads$ctrl_v_reac) %>% arrange(desc(abs(PC1)))

# calculate 'optimum' number of PCs to use via horn and elbow methods, plot scree
# also how many components for 80%
optimumScree <- function(pca, svcm){
  horn <- parallelPCA(svcm)
  elbow <- findElbowPoint(pca$variance)
  comp <- which(cumsum(pca$variance) > 80)[1]
  plt <- screeplot(pca,
            components = getComponents(pca, 1:8),
            vline = c(comp, horn$n, elbow), 
            axisLabSize = 8) +
    geom_label(aes(x = horn$n, y = 50,
                   label = 'Horn\'s', vjust = -3, size = 8)) +
    geom_label(aes(x = elbow, y = 50,
                   label = 'Elbow method', vjust = -3, size = 8))  +
    geom_label(aes(x = comp, y = 50,
                   label = '80% explained', vjust = -3, size = 8)) 
  return(plt)
}

for (name in names(svcms)){
  optimumScree(pcas[[name]], svcms[[name]])
}  

# number of components that need to be included to explain X% variance
#ncomponents <- lapply(pcas, function(p) which(cumsum(p$variance) > 80)[1])

