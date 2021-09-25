# here we import the TPS file of landmarks and do geometric morphometry
# https://www-nature-com.offcampus.lib.washington.edu/articles/s41596-020-0347-z

library(tidyverse)
library(geomorph)
library(Morpho)
library(mclust)
library(cluster)
library(factoextra)

lands <- readland.tps(
  here::here("data/landmarks.tps"),
  specID = "ID")

# drop outliers  identified in the previous script
outlier_ids <- c(1152,  606, 1147, 1125,  756 ,1023,  457,
                 813,  415,  877,  677, 1055,  239,  811,  601,  861,
                  794,  103 , 598 ,1110 , 906 ,  59 ,   4)

lands_no_outliers <- lands[, , -outlier_ids]
lands_no_outliers_names <- dimnames(lands_no_outliers)[[3]]

sup.geo <- gpagen(lands_no_outliers,  ProcD=F, Proj = T)

# Principal component analysis (PCA) using a base R function
pca.geo <- prcomp(two.d.array(sup.geo$coord))

# TODO: find a meaningful way to break this up into groups
# pass these groups onto the PCA for drawing ellipses

# this block of code below shows one way to discretize cortex percentage
# into three groups. We need to explore other cutoffs, and other grouping
# variables besides dorsal cortex. This code below is a good template for
# starting these explorations
ma_flakes_tr_with_groups <-
  ma_flakes_tr %>%
  slice(-outlier_ids) %>%
  mutate(dorsal_cortext_group =
           case_when(
             between(dorsal_cortex_perc, 0, 10) ~ "group 1",
             between(dorsal_cortex_perc, 11, 60) ~ "group 2",
             between(dorsal_cortex_perc, 61, 100) ~ "group 3"
           ))

# ggplot2-based plot of the PCA results using a factoextra function
fviz_pca_ind(pca.geo,
             geom.ind = "point", # show points without labels
            #  line below sets the grouping column for colouring the points
             habillage = ma_flakes_tr_with_groups$dorsal_cortext_group,
             palette = "npg", #color palette, adjust to preferences
             addEllipses = T, #adds ellipses or convex hulls
             ellipse.type = "convex", #convex hulls
             mean.point = F, #don’t show average point
             legend.title = "SP") +  #adjust accordingly)
             coord_fixed(ratio = 1) #adjustable for variation-weighted plots

# let's look at all the PCs and how much each original # variable contributes to them

#view how much variation is described by each PC (ordered from PC1 to last PC)
enframe(100*pca.geo$sdev^2/sum(pca.geo$sdev^2))

#get meaningful PCs
getMeaningfulPCs(pca.geo$sdev^2, n = length(lands_no_outliers_names))

#estimate contribution (%) of landmarks on PC1
PC1_contrib <- fviz_contrib(pca.geo, choice = "var", axes = 1)
# take a look at the plot
PC1_contrib

LM_contrib_PC1 <- data.frame(LM = rep(NA, nrow(PC1_contrib$data)/2),
                             PC1_contrib = rep(NA, nrow(PC1_contrib$data)/2))
for (i in 1:(nrow(PC1_contrib$data)/2)) {
  LM_contrib_PC1$PC1_contrib[i] <- sum(PC1_contrib$data[(i*2-1):(i*2),]$contrib)
  LM_contrib_PC1$LM[i] <- as.character(i)
}

# TODO: decode the x-axis labels on the scree plot so we can describe in
# terms of actual flake dimension, write some sentences to describe and interpret

# estimate contribution (%) of landmarks on PC2
PC2_contrib <- fviz_contrib(pca.geo, choice = "var", axes = 2)
# take a look at the plot
PC2_contrib

LM_contrib_PC2 <- data.frame(LM = rep(NA, nrow(PC2_contrib$data)/2),
                             PC2_contrib = rep(NA, nrow(PC2_contrib$data)/2))
for (i in 1:(nrow(PC2_contrib$data)/2)) {
  LM_contrib_PC2$PC2_contrib[i] <- sum(PC2_contrib$data[(i*2-1):(i*2),]$contrib)
  LM_contrib_PC2$LM[i] <- as.character(i)
}

# examine landmarks with major contribution (%) to PC1 and PC2
important_LMs_PC1 <- LM_contrib_PC1[LM_contrib_PC1$PC1_contrib > 100/(nrow(PC1_contrib$data)/2),]
important_LMs_PC1[order(decreasing = T, important_LMs_PC1$PC1_contrib),]
important_LMs_PC2 <- LM_contrib_PC2[LM_contrib_PC2$PC2_contrib > 100/(nrow(PC2_contrib$data)/2),]
# here is a table that summarises the most important landmarks in the PCA
important_LMs_PC2[order(decreasing = T, important_LMs_PC2$PC2_contrib),]

