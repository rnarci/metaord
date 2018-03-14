########################## Stocker les clusterings pour chaque matrice, nombre de clusters et nombre 
########################## de plus proches voisins

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


source('~/metaord/utils.R')

size_fraction3 = "0.8-5"

import_data(size_fraction3, samples = NULL)

n = dim(jaccard_abundance)[1]

kNN1 = matrix(NA, nrow = n, ncol = n)
kNN2 = matrix(NA, nrow = n, ncol = n)
kNN3 = matrix(NA, nrow = n, ncol = n)
kNN4 = matrix(NA, nrow = n, ncol = n)
kNN5 = matrix(NA, nrow = n, ncol = n)
kNN6 = matrix(NA, nrow = n, ncol = n)
kNN7 = matrix(NA, nrow = n, ncol = n)
kNN8 = matrix(NA, nrow = n, ncol = n)
kNN9 = matrix(NA, nrow = n, ncol = n)
kNN10 = matrix(NA, nrow = n, ncol = n)
kNN11 = matrix(NA, nrow = n, ncol = n)
kNN12 = matrix(NA, nrow = n, ncol = n)

res1 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res2 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res3 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res4 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res5 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res6 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res7 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res8 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res9 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res10 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res11 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res12 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)

kNN = list(kNN1,kNN2,kNN3,kNN4,kNN5,kNN6,kNN7,kNN8,kNN9,kNN10,kNN11,kNN12)
res = list(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12)
matrices.list = c("jaccard_abundance", "ab_jaccard_abundance", "braycurtis_abundance", 
  "ab_ochiai_abundance", "ab_sorensen_abundance", "simka_jaccard_abundance", 
  "chord_prevalence", "jaccard_prevalence", "kulczynski_prevalence", 
  "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence")

a = 0

for(k in 1:(round(n/2))){
  for(num_matrix in 1:12){
    kNN[[num_matrix]] <- make.kNNG(get(matrices.list[num_matrix]), k = k, symm = TRUE, weight = FALSE)
    for(l in 2:15){
      res[[num_matrix]][14*(k-1)+l-1,1] = k 
      res[[num_matrix]][14*(k-1)+l-1,2] = l
      res[[num_matrix]][14*(k-1)+l-1,3:(n+2)] = spectral.clustering.new(kNN[[num_matrix]], normalised = TRUE, score = FALSE, K = l, adj = FALSE)
    }
  }
  a = a + 14
  print(a)
}

for(num_matrix in 1:12){
  res[[num_matrix]] = as.data.frame(res[[num_matrix]])
}

write.table(x = res[[1]], file = "Mat1_spectral")
write.table(x = res[[2]], file = "Mat2_spectral")
write.table(x = res[[3]], file = "Mat3_spectral")
write.table(x = res[[4]], file = "Mat4_spectral")
write.table(x = res[[5]], file = "Mat5_spectral")
write.table(x = res[[6]], file = "Mat6_spectral")
write.table(x = res[[7]], file = "Mat7_spectral")
write.table(x = res[[8]], file = "Mat8_spectral")
write.table(x = res[[9]], file = "Mat9_spectral")
write.table(x = res[[10]], file = "Mat10_spectral")
write.table(x = res[[11]], file = "Mat11_spectral")
write.table(x = res[[12]], file = "Mat12_spectral")





















########################## Stocker les clusterings pour chaque matrice, nombre de clusters et dimension 
########################## de l'espace dans lequel on projette les echantillons


rm(list=objects())

library(fcd)
library(loe)
library(clues)
library(cccd)
library(fossil)
library(igraph)


source('~/metaord/utils.R')

size_fraction3 = "0.8-5"

import_data(size_fraction3, samples = NULL)

n = dim(jaccard_abundance)[1]

res1 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res2 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res3 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res4 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res5 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res6 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res7 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res8 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res9 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res10 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res11 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
res12 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)

res = list(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12)

a = 0

for(k in 1:(round(n/2))){
    fit1 = cmdscale(jaccard_abundance, k=k)
    fit2 = cmdscale(ab_jaccard_abundance, k=k)
    fit3 = cmdscale(braycurtis_abundance, k=k)
    fit4 = cmdscale(ab_ochiai_abundance, k=k)
    fit5 = cmdscale(ab_sorensen_abundance, k=k)
    fit6 = cmdscale(simka_jaccard_abundance, k=k)
    fit7 = cmdscale(chord_prevalence, k=k)
    fit8 = cmdscale(jaccard_prevalence, k=k)
    fit9 = cmdscale(kulczynski_prevalence, k=k)
    fit10 = cmdscale(ochiai_prevalence, k=k)
    fit11 = cmdscale(whittaker_prevalence, k=k)
    fit12 = cmdscale(simka_jaccard_prevalence, k=k)
    for(l in 2:15){
      for(num_matrix in 1:12){
        res[[num_matrix]][14*(k-1)+l-1,1] = k 
        res[[num_matrix]][14*(k-1)+l-1,2] = l
        res[[num_matrix]][14*(k-1)+l-1,3:(n+2)] = kmeans(x = get(paste("fit",num_matrix,sep="")), centers = l, nstart = 1000, iter.max = 20)$cluster
      }
    }
  a = a + 14
  print(a)
}

for(num_matrix in 1:12){
  res[[num_matrix]] = as.data.frame(res[[num_matrix]])
}

write.table(x = res[[1]], file = "Mat1_kmeans")
write.table(x = res[[2]], file = "Mat2_kmeans")
write.table(x = res[[3]], file = "Mat3_kmeans")
write.table(x = res[[4]], file = "Mat4_kmeans")
write.table(x = res[[5]], file = "Mat5_kmeans")
write.table(x = res[[6]], file = "Mat6_kmeans")
write.table(x = res[[7]], file = "Mat7_kmeans")
write.table(x = res[[8]], file = "Mat8_kmeans")
write.table(x = res[[9]], file = "Mat9_kmeans")
write.table(x = res[[10]], file = "Mat10_kmeans")
write.table(x = res[[11]], file = "Mat11_kmeans")
write.table(x = res[[12]], file = "Mat12_kmeans")















########################### Etudes de la stabilite a travers le calcul des RI


########### Spectral clustering 

rm(list=objects())

# library(fcd)
# library(loe)
# library(clues)
# library(cccd)
# library(fossil)
# library(igraph)


source('~/metaord/utils.R')

size_fraction3 = "0.8-5"

import_data(size_fraction3, samples = NULL)

n = dim(jaccard_abundance)[1]

mat1 = read.table(file = "Mat1_spectral")
mat2 = read.table(file = "Mat2_spectral")
mat3 = read.table(file = "Mat3_spectral")
mat4 = read.table(file = "Mat4_spectral")
mat5 = read.table(file = "Mat5_spectral")
mat6 = read.table(file = "Mat6_spectral")
mat7 = read.table(file = "Mat7_spectral")
mat8 = read.table(file = "Mat8_spectral")
mat9 = read.table(file = "Mat9_spectral")
mat10 = read.table(file = "Mat10_spectral")
mat11 = read.table(file = "Mat11_spectral")
mat12 = read.table(file = "Mat12_spectral")

mat = list(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8,mat9,mat10,mat11,mat12)

for(num_matrix in 1:12){
  mat[[num_matrix]] = as.matrix(mat[[num_matrix]])
}


RI1 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI2 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI3 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI4 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI5 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI6 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI7 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI8 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI9 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI10 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI11 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI12 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)

RI = list(RI1,RI2,RI3,RI4,RI5,RI6,RI7,RI8,RI9,RI10,RI11,RI12)

for(num_matrix in 1:12){
  a = 1
  for(k1 in 1:(round(n/2)-1)){
    for(k2 in (k1+1):round(n/2)){
      for(l in 2:15){
        RI[[num_matrix]][a,1] = l
        RI[[num_matrix]][a,2] = k1
        RI[[num_matrix]][a,3] = k2
        RI[[num_matrix]][a,4] = rand.index(mat[[num_matrix]][14*(k1-1)+l-1,3:(n+2)], mat[[num_matrix]][14*(k2-1)+l-1,3:(n+2)])
        a = a + 1
        print(a)
      }
    }
  }
}


for(num_matrix in 1:12){
  RI[[num_matrix]] = as.data.frame(RI[[num_matrix]])
}
 
# write.table(x = RI[[1]], file = "RI1_spectral")
# write.table(x = RI[[2]], file = "RI2_spectral")
# write.table(x = RI[[3]], file = "RI3_spectral")
# write.table(x = RI[[4]], file = "RI4_spectral")
# write.table(x = RI[[5]], file = "RI5_spectral")
# write.table(x = RI[[6]], file = "RI6_spectral")
# write.table(x = RI[[7]], file = "RI7_spectral")
# write.table(x = RI[[8]], file = "RI8_spectral")
# write.table(x = RI[[9]], file = "RI9_spectral")
# write.table(x = RI[[10]], file = "RI10_spectral")
# write.table(x = RI[[11]], file = "RI11_spectral")
# write.table(x = RI[[12]], file = "RI12_spectral")

# res <- lapply(1:12, function(i) { 
#   res <- read.table(paste0("RI", i, "_spectral"))
#   colnames(res) <- c("Nb_Cluster", "Diff.k", "RI")
#   res %>% mutate(Matrice.Nb = i)
#   }) %>% 
#   bind_rows() %>% 
#   mutate(Distance = matrices.list[Matrice.Nb])

matrices.list = c("jaccard_abundance", "ab_jaccard_abundance", "braycurtis_abundance", 
                  "ab_ochiai_abundance", "ab_sorensen_abundance", "simka_jaccard_abundance", 
                  "chord_prevalence", "jaccard_prevalence", "kulczynski_prevalence", 
                  "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence")

res <- lapply(1:12, function(i) { 
  res <- RI[[i]]
  colnames(res) <- c("Nb_Cluster", "k1", "k2", "RI")
  res %>% mutate(Matrice.Nb = i)
}) %>% 
  bind_rows() %>% 
  mutate(Distance = matrices.list[Matrice.Nb])

fwrite(res, file = "RI_spectral_all")

########### kmeans

rm(list=objects())

source('~/metaord/utils.R')

size_fraction3 = "0.8-5"

import_data(size_fraction3, samples = NULL)

n = dim(jaccard_abundance)[1]

mat1 = read.table(file = "Mat1_kmeans")
mat2 = read.table(file = "Mat2_kmeans")
mat3 = read.table(file = "Mat3_kmeans")
mat4 = read.table(file = "Mat4_kmeans")
mat5 = read.table(file = "Mat5_kmeans")
mat6 = read.table(file = "Mat6_kmeans")
mat7 = read.table(file = "Mat7_kmeans")
mat8 = read.table(file = "Mat8_kmeans")
mat9 = read.table(file = "Mat9_kmeans")
mat10 = read.table(file = "Mat10_kmeans")
mat11 = read.table(file = "Mat11_kmeans")
mat12 = read.table(file = "Mat12_kmeans")

mat = list(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8,mat9,mat10,mat11,mat12)

for(num_matrix in 1:12){
  mat[[num_matrix]] = as.matrix(mat[[num_matrix]])
}

RI1 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI2 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI3 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI4 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI5 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI6 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI7 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI8 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI9 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI10 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI11 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
RI12 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)

RI = list(RI1,RI2,RI3,RI4,RI5,RI6,RI7,RI8,RI9,RI10,RI11,RI12)

for(num_matrix in 1:12){
  a = 1
  for(k1 in 1:(round(n/2)-1)){
    for(k2 in (k1+1):round(n/2)){
      for(l in 2:15){
        RI[[num_matrix]][a,1] = l
        RI[[num_matrix]][a,2] = k1
        RI[[num_matrix]][a,3] = k2
        RI[[num_matrix]][a,4] = rand.index(mat[[num_matrix]][14*(k1-1)+l-1,3:(n+2)], mat[[num_matrix]][14*(k2-1)+l-1,3:(n+2)])
        a = a + 1
        print(a)
      }
    }
  }
}


for(num_matrix in 1:12){
  RI[[num_matrix]] = as.data.frame(RI[[num_matrix]])
}

matrices.list = c("jaccard_abundance", "ab_jaccard_abundance", "braycurtis_abundance", 
                  "ab_ochiai_abundance", "ab_sorensen_abundance", "simka_jaccard_abundance", 
                  "chord_prevalence", "jaccard_prevalence", "kulczynski_prevalence", 
                  "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence")

res <- lapply(1:12, function(i) { 
  res <- RI[[i]]
  colnames(res) <- c("Nb_Cluster", "d1", "d2", "RI")
  res %>% mutate(Matrice.Nb = i)
}) %>% 
  bind_rows() %>% 
  mutate(Distance = matrices.list[Matrice.Nb])

fwrite(res, file = "RI_kmeans_all")

# write.table(x = RI[[1]], file = "RI1_kmeans")
# write.table(x = RI[[2]], file = "RI2_kmeans")
# write.table(x = RI[[3]], file = "RI3_kmeans")
# write.table(x = RI[[4]], file = "RI4_kmeans")
# write.table(x = RI[[5]], file = "RI5_kmeans")
# write.table(x = RI[[6]], file = "RI6_kmeans")
# write.table(x = RI[[7]], file = "RI7_kmeans")
# write.table(x = RI[[8]], file = "RI8_kmeans")
# write.table(x = RI[[9]], file = "RI9_kmeans")
# write.table(x = RI[[10]], file = "RI10_kmeans")
# write.table(x = RI[[11]], file = "RI11_kmeans")
# write.table(x = RI[[12]], file = "RI12_kmeans")


########################################## Representation des resultats

setwd(dir="~/metaord/RI_results")

mean.plus.sd <- function(x){
  ## return(min(mean(x) + sd(x), 1))
  quantile(x, 0.9)
}

mean.moins.sd <- function(x){
  ## return(max(mean(x) - sd(x), 0))
  quantile(x, 0.1)
}




############### Spectral clustering

res = fread("RI_spectral_all")
res$Diff.k = abs(res$k1 - res$k2)

# ggplot(res %>% filter(Distance == "jaccard_abundance"), aes(x = Diff.k, y = RI)) +
#   geom_point(alpha = 0.3, size = 0.5) +
#   # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
#   facet_wrap(~Nb_Cluster)

ggplot(res, aes(x = Diff.k, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)

ggplot(res %>% filter(k1 > 10 & k2 > 10), aes(x = Diff.k, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)

ggplot(res %>% filter(k1 > 30 & k2 > 30), aes(x = Diff.k, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)

ggplot(res %>% filter(Diff.k <= 20), aes(x = Diff.k, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)

# ggplot(res %>% filter(Distance == "jaccard_abundance"), aes(x = Diff.k, y = RI)) +
#   stat_summary(fun.y = "mean", colour = "red", size = 1, geom = "path") +
#   facet_wrap(~Nb_Cluster)


# ggplot(res, aes(x = Diff.k, y = RI, color = Distance)) +
#   stat_summary(fun.y = "mean", size = 1, geom = "path") +
#   stat_summary(fun.y = "mean.moins.sd", size = 0.5, geom = "path", colour = "black") +
#   stat_summary(fun.y = "mean.plus.sd", size = 0.5, geom = "path", colour = "black") +
#   facet_grid(Distance~Nb_Cluster)
# 
# ggplot(res %>% filter(k1 > 10 & k2 > 10), aes(x = Diff.k, y = RI, color = Distance)) +
#   stat_summary(fun.y = "mean", size = 1, geom = "path") +
#   stat_summary(fun.y = "mean.moins.sd", size = 0.5, geom = "path", colour = "black") +
#   stat_summary(fun.y = "mean.plus.sd", size = 0.5, geom = "path", colour = "black") +
#   facet_grid(Distance~Nb_Cluster)
# 
# ggplot(res %>% filter(k1 > 30 & k2 > 30), aes(x = Diff.k, y = RI, color = Distance)) +
#   stat_summary(fun.y = "mean", size = 1, geom = "path") +
#   stat_summary(fun.y = "mean.moins.sd", size = 0.5, geom = "path", colour = "black") +
#   stat_summary(fun.y = "mean.plus.sd", size = 0.5, geom = "path", colour = "black") +
#   facet_grid(Distance~Nb_Cluster)
# 
# p <- ggplot(res %>% filter(Diff.k <= 20), aes(x = Diff.k, y = RI, color = Distance)) +
#   stat_summary(fun.y = "mean", size = 1, geom = "path") +
#   stat_summary(fun.y = "mean.moins.sd", size = 0.5, geom = "path", colour = "black") +
#   stat_summary(fun.y = "mean.plus.sd", size = 0.5, geom = "path", colour = "black") +
#   facet_grid(Distance~Nb_Cluster)
# plot(p)
# ggsave(p, file = "myfile.png", width = 15, height = 10, dpi = 100)

############## kmeans

res = fread("RI_kmeans_all")

res$Diff.d = abs(res$d1 - res$d2)

# ggplot(res %>% filter(Distance == "jaccard_abundance"), aes(x = Diff.d, y = RI)) +
#   geom_point(alpha = 0.3, size = 0.5) +
#   # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
#   facet_wrap(~Nb_Cluster)

ggplot(res, aes(x = Diff.d, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)

ggplot(res %>% filter(d1 > 10 & d2 > 10), aes(x = Diff.d, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)

ggplot(res %>% filter(d1 > 30 & d2 > 30), aes(x = Diff.d, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)

ggplot(res %>% filter(Diff.d <= 20), aes(x = Diff.d, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)

# ggplot(res %>% filter(Distance == "jaccard_abundance"), aes(x = Diff.d, y = RI)) +
#   stat_summary(fun.y = "mean", colour = "red", size = 1, geom = "path") +
#   facet_wrap(~Nb_Cluster)

# ggplot(res, aes(x = Diff.d, y = RI, color = Distance)) +
#   stat_summary(fun.y = "mean", size = 1, geom = "path") +
#   stat_summary(fun.y = "mean.moins.sd", size = 0.5, geom = "path", colour = "black") +
#   stat_summary(fun.y = "mean.plus.sd", size = 0.5, geom = "path", colour = "black") +
#   facet_grid(Distance~Nb_Cluster)
# 
# ggplot(res %>% filter(d1 > 10 & d2 > 10), aes(x = Diff.d, y = RI, color = Distance)) +
#   stat_summary(fun.y = "mean", size = 1, geom = "path") +
#   stat_summary(fun.y = "mean.moins.sd", size = 0.5, geom = "path", colour = "black") +
#   stat_summary(fun.y = "mean.plus.sd", size = 0.5, geom = "path", colour = "black") +
#   facet_grid(Distance~Nb_Cluster)
# 
# ggplot(res %>% filter(d1 > 30 & d2 > 30), aes(x = Diff.d, y = RI, color = Distance)) +
#   stat_summary(fun.y = "mean", size = 1, geom = "path") +
#   stat_summary(fun.y = "mean.moins.sd", size = 0.5, geom = "path", colour = "black") +
#   stat_summary(fun.y = "mean.plus.sd", size = 0.5, geom = "path", colour = "black") +
#   facet_grid(Distance~Nb_Cluster)
# 
# ggplot(res %>% filter(Diff.d <= 20), aes(x = Diff.d, y = RI, color = Distance)) +
#   stat_summary(fun.y = "mean", size = 1, geom = "path") +
#   stat_summary(fun.y = "mean.moins.sd", size = 0.5, geom = "path", colour = "black") +
#   stat_summary(fun.y = "mean.plus.sd", size = 0.5, geom = "path", colour = "black") +
#   facet_grid(Distance~Nb_Cluster)

############################ Superposer moyenne des RI (kmeans et spectral)

res1 = fread("RI_kmeans_all") %>% rename(param1 = d1, param2 = d2)
res2 = fread("RI_spectral_all") %>% rename(param1 = k1, param2 = k2)
#res1$Diff = abs(res1$d1 - res1$d2)
#res2$Diff = abs(res2$k1 - res2$k2)
res <- bind_rows(res1 %>% mutate(Type = "K-means"), 
                 res2 %>% mutate(Type = "Spectral")) %>% 
  mutate(Diff = abs(param1 - param2))

res1 <- res %>% group_by(Nb_Cluster, Distance, Diff, Type) %>% 
  summarize(mean.RI   = mean(RI), 
            sd.RI     = sd(RI), 
            q90       = quantile(RI, 0.9), 
            median.RI = quantile(RI, 0.5), 
            q10       = quantile(RI, 0.1)) %>% 
  gather(key = "Summary", value = "RI", mean.RI:q10)

p = ggplot(res1 %>% filter(Summary %in% c("q10", "median.RI", "q90")), 
       aes(x = Diff, y = RI, group = interaction(Summary, Type), color = Type)) + 
  geom_line(aes(linetype = (Summary != "median.RI")), size = 0.5) +
  scale_linetype_discrete(guide = FALSE) + 
  facet_grid(Distance~Nb_Cluster) 
ggsave(p, file = "all_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/summary_RI", width = 15, height = 10, dpi = 100)

res1 <- res %>% filter(param1 > 10 & param2 > 10) %>% group_by(Nb_Cluster, Distance, Diff, Type) %>% 
  summarize(mean.RI   = mean(RI), 
            sd.RI     = sd(RI), 
            q90       = quantile(RI, 0.9), 
            median.RI = quantile(RI, 0.5), 
            q10       = quantile(RI, 0.1)) %>% 
  gather(key = "Summary", value = "RI", mean.RI:q10)

p = ggplot(res1 %>% filter(Summary %in% c("q10", "median.RI", "q90")), 
           aes(x = Diff, y = RI, group = interaction(Summary, Type), color = Type)) + 
  geom_line(aes(linetype = (Summary != "median.RI")), size = 0.5) +
  scale_linetype_discrete(guide = FALSE) + 
  facet_grid(Distance~Nb_Cluster) 
ggsave(p, file = "seuil1_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/summary_RI", width = 15, height = 10, dpi = 100)

res1 <- res %>% filter(param1 > 30 & param2 > 30) %>% group_by(Nb_Cluster, Distance, Diff, Type) %>% 
  summarize(mean.RI   = mean(RI), 
            sd.RI     = sd(RI), 
            q90       = quantile(RI, 0.9), 
            median.RI = quantile(RI, 0.5), 
            q10       = quantile(RI, 0.1)) %>% 
  gather(key = "Summary", value = "RI", mean.RI:q10)

p = ggplot(res1 %>% filter(Summary %in% c("q10", "median.RI", "q90")), 
           aes(x = Diff, y = RI, group = interaction(Summary, Type), color = Type)) + 
  geom_line(aes(linetype = (Summary != "median.RI")), size = 0.5) +
  scale_linetype_discrete(guide = FALSE) + 
  facet_grid(Distance~Nb_Cluster) 
ggsave(p, file = "seuil2_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/summary_RI", width = 15, height = 10, dpi = 100)

res1 <- res %>% filter(Diff <= 20) %>% group_by(Nb_Cluster, Distance, Diff, Type) %>% 
  summarize(mean.RI   = mean(RI), 
            sd.RI     = sd(RI), 
            q90       = quantile(RI, 0.9), 
            median.RI = quantile(RI, 0.5), 
            q10       = quantile(RI, 0.1)) %>% 
  gather(key = "Summary", value = "RI", mean.RI:q10)

p = ggplot(res1 %>% filter(Summary %in% c("q10", "median.RI", "q90")), 
           aes(x = Diff, y = RI, group = interaction(Summary, Type), color = Type)) + 
  geom_line(aes(linetype = (Summary != "median.RI")), size = 0.5) +
  scale_linetype_discrete(guide = FALSE) + 
  facet_grid(Distance~Nb_Cluster) 
ggsave(p, file = "seuil3_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/summary_RI", width = 15, height = 10, dpi = 100)


# ggplot() + 
#   stat_summary(data = res1 , aes(x = Diff, y = RI), fun.y = "mean", size = 0.5, geom = "path", colour = "red") +
#   stat_summary(data = res2, aes(x = Diff, y = RI), fun.y = "mean", size = 0.5, geom = "path", colour = "blue") +
#   facet_grid(Distance~Nb_Cluster) 
# 
# ggplot() + 
#   stat_summary(data = res1 %>% filter(d1 > 10 & d2 > 10), aes(x = Diff, y = RI), fun.y = "mean", size = 0.5, geom = "path", colour = "red") +
#   stat_summary(data = res2 %>% filter(k1 > 10 & k2 > 10), aes(x = Diff, y = RI), fun.y = "mean", size = 0.5, geom = "path", colour = "blue") +
#   facet_grid(Distance~Nb_Cluster) 
# 
# ggplot() + 
#   stat_summary(data = res1 %>% filter(d1 > 30 & d2 > 30), aes(x = Diff, y = RI), fun.y = "mean", size = 0.5, geom = "path", colour = "red") +
#   stat_summary(data = res2 %>% filter(k1 > 30 & k2 > 30), aes(x = Diff, y = RI), fun.y = "mean", size = 0.5, geom = "path", colour = "blue") +
#   facet_grid(Distance~Nb_Cluster) 
# 
# ggplot() + 
#   stat_summary(data = res1 %>% filter(Diff <= 20), aes(x = Diff, y = RI), fun.y = "mean", size = 0.5, geom = "path", colour = "red") +
#   stat_summary(data = res2 %>% filter(Diff <= 20), aes(x = Diff, y = RI), fun.y = "mean", size = 0.5, geom = "path", colour = "blue") +
#   facet_grid(Distance~Nb_Cluster) 
