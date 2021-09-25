# here we import the TPS file of landmarks and do geometric morphometry
# https://www-nature-com.offcampus.lib.washington.edu/articles/s41596-020-0347-z

# This script is for initial exploration and outlier detection. The
# next script file is the one to look at for results to describe and
# interpret

library(tidyverse)
library(geomorph)
library(Morpho)
library(mclust)
library(cluster)
library(factoextra)

lands <- readland.tps(
  here::here("data/landmarks.tps"),
  specID = "ID")

lands.names <- dimnames(lands)[[3]]

sup.geo <- gpagen(lands,  ProcD=F, Proj = T)

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

# consider removing these and repeating the GPA/PCA
outliers <- plotOutliers(sup.geo$coords) #generate an outlier plot
# get outlier IDs
outlier_ids <- unname(outliers)[1:(length(outliers)*0.02)]


#if a potential outlier is found, compare it to the mean shape of the data set
shape1 <- sup.geo$consensus #mean shape
shape2 <- sup.geo$coords[,,1]
#adjust the value so that it fits the number that corresponds to the outlier
#generate deformation grid of mean shape and potential outlier to check for errors
plotRefToTarget(shape1, shape2, method = "TPS", mag = 2)

#----------------------------------------
