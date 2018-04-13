rm(list=objects())

library(fcd)
library(loe)
library(clues)
library(cccd)
library(fossil)
library(igraph)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gdata)
library(combinat)

source('~/metaord/utils.R')

############################################# Import Genoscope data

oldwd <- getwd() 

setwd(dir="~/Bureau/CDDrnarci/Donnees")
data1 = read.xls(xls="Genocenoses_env_parameters_all_tara.xlsx", method = "csv")
# data2 = read.xls(xls="Genocenoses_env_parameters_all_woa.xlsx", method = "csv") -> same clustering
data.0_0.2 = data1[which((as.character(data1$Fraction) == "0-0.2") == TRUE),c(1,3)] 
data.0.22_3 = data1[which((as.character(data1$Fraction) == "0.22-3") == TRUE),c(1,3)] 
data.0.8_5 = data1[which((as.character(data1$Fraction) == "0.8-5") == TRUE),c(1,3)] 
data.5_20 = data1[which((as.character(data1$Fraction) == "May-20") == TRUE),c(1,3)] 
data.20_180 = data1[which((as.character(data1$Fraction) == "20-180") == TRUE),c(1,3)] 
data.180_2000 = data1[which((as.character(data1$Fraction) == "180-2000") == TRUE),c(1,3)] 

rownames(data.0_0.2) = data.0_0.2$Station
rownames(data.0.22_3) = data.0.22_3$Station
rownames(data.0.8_5) = data.0.8_5$Station
rownames(data.5_20) = data.5_20$Station
rownames(data.20_180) = data.20_180$Station
rownames(data.180_2000) = data.180_2000$Station

setwd(oldwd)

# size_fraction = "0.8-5"
# import_data(size_fraction, samples = rownames(data.0.8_5))
# n = dim(jaccard_abundance)[1]

########################## To stock the clusterings for each matrix, the number of clusters and either
########################## the dimension of the space in which we project the samples or the number 
########################## of nearest neighbours 


# get_data_RI()


######################################################## Plots of the results

# type = "spectral" # spectral clustering based on ordinal distances or kmeans on the new space of 
# # the projected points
# method = "RI_plot" # plot of all of the RIs or a summary (median, q10, q90)
# threshold = c(11,n-1) # to threshold the parameters (k (kNNG) or d (cmdscale))
#   
# p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
# plot(p)

# ggsave(p, file = "seuil3_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/summary_RI", width = 15, height = 10, dpi = 100)

######################################################## Comparison with the clusterings obtained by Genoscope

comparison = compare_clusterings()

size_fraction = "0-0.2"
frac = which(comparison$Fraction == size_fraction)

p = ggplot(comparison[frac,], aes(x = Param, y = RI, color = Type)) +
  geom_point(alpha = 1, size = 1) +
  scale_y_continuous(limits = c(0,1)) +
  # ggtitle(size_fraction) +
  facet_wrap(~Distance)
plot(p)

p = ggplot(comparison, aes(x = Distance, y = RI, color = Type)) +
  geom_boxplot() +
  # scale_y_continuous(limits = c(0,1)) +
  # ggtitle(size_fraction) +
  facet_wrap(~Fraction)
plot(p)

######################################################## Identify the samples which are clustered in different group :
######################################################## K-means clustering vs Genoscope clustering

comparison = compare_clusterings()

size_fraction = "0.8-5"
id_samples = identify_samples_ex(size_fraction)
# id_samples = identify_samples(size_fraction)
# Kmeans_clustering = id_samples$Kmeans_clustering
# Kmeans_initial = id_samples$Kmeans_initial
# Genoscope_clustering = id_samples$Genoscope_clustering
all_permutations = id_samples$all_permutations
score = id_samples$score
max(score)
which(score==max(score))
all_permutations[472,]
all_permutations[489,]

score
rand.index(Kmeans_clustering,Kmeans_initial)
rand.index(Kmeans_clustering,Genoscope_clustering)
Kmeans_clustering
Genoscope_clustering

which(Kmeans_clustering!=Genoscope_clustering)

############################################################# Eigenvalues study

size_fraction = "180-2000"
eigenvalues_study(size_fraction)
