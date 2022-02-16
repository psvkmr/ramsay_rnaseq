library(PCAtools)
library(DESeq2)
library(tidyverse)
library(vsn)

#https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html

setwd('C:/Users/Prasanth/Documents/ramsey/')

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

# in this analysis, only the 'NAT' samples are going to be used, 
# so both the md and cm will be subsetted for these
smpls <- colnames(zcm)[grep('NAT', colnames(zcm))]
nat.cm <- zcm[, smpls]
nat.md <- md[smpls, ]


# sample groups
smpl.grps <- list(cll_ln = rownames(nat.md[nat.md$sample_group == 'CLL LN\n', ]), 
                  cll_pb = rownames(nat.md[nat.md$sample_group == 'CLL PB\n', ]), 
                  ctrl = rownames(nat.md[nat.md$sample_group == 'Control\n', ]), 
                  reac = rownames(nat.md[nat.md$sample_group == 'Reactive LN', ]))

# cms subset by sample group comparisons
zcms <- list(all = zcm,
             ctrl_v_cll_ln = zcm[, c(smpl.grps$ctrl, smpl.grps$cll_ln)], 
             ctrl_v_cll_pb = zcm[, c(smpl.grps$ctrl, smpl.grps$cll_pb)], 
             ctrl_v_reac = zcm[, c(smpl.grps$ctrl, smpl.grps$reac)])
#lapply(zcms, dim)

# cms filtered for low counts specifically in comparison subset
zcms <- lapply(zcms, function(mat) mat[rowSums(mat) > (((dim(mat)[2])/2)-1), ])
#lapply(zcms, dim)

# normalised subset cms
ncms <- list(all = ncm,
             ctrl_v_cll_ln = ncm[, c(smpl.grps$ctrl, smpl.grps$cll_ln)], 
             ctrl_v_cll_pb = ncm[, c(smpl.grps$ctrl, smpl.grps$cll_pb)], 
             ctrl_v_reac = ncm[, c(smpl.grps$ctrl, smpl.grps$reac)])
#lapply(ncms, dim)

