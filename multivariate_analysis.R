rm(list=objects())

library(fcd)
library(loe)
library(clues)
library(cccd)
library(fossil)
library(igraph)

############################################################ Import data

size_fraction = "20-180"
import_data(size_fraction)

############################## Design

design = read.table(file="param_bioadvection.csv",sep="",header=TRUE)

library(vegan)

adonis(jaccard_abundance[1:130,] ~ design[1:130,2], data=design, permutations = 999)
