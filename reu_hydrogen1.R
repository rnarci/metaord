rm(list=objects())

library(fcd)
library(loe)
library(clues)
library(cccd)
library(fossil)
library(igraph)

############################################################ Source custom scripts

source('~/metaord/utils.R')

############################################################ Import data


size_fraction = "0.8-5"
import_data(size_fraction, samples = NULL)

############################################################ Clustering

############################## MDS + k-means

###### Abundance matrices

# l = 5
# iter = 1
# RI = matrix(NA,9,3)
# 
# for(k1 in c(4,15,40))
# {
#   fit1 = cmdscale(jaccard_abundance, eig=TRUE, k=k1)
#   fit2 = cmdscale(ab_sorensen_abundance, eig=TRUE, k=k1)
#   for(k2 in c(4,15,40)){
#     fitbis1 = cmdscale(jaccard_abundance, eig=TRUE, k=k2)
#     fitbis2 = cmdscale(ab_sorensen_abundance,eig=TRUE, k=k2)
# 
#     res1 = kmeans(fit1$points,l,nstart = 25)$cluster
#     res2 = kmeans(fit2$points,l,nstart = 25)$cluster
#     resbis1 = kmeans(fitbis1$points,l,nstart = 25)$cluster
#     resbis2 = kmeans(fitbis2$points,l,nstart = 25)$cluster
# 
#     RI[iter,1] = rand.index(get(paste("res",1,sep="")),get(paste("resbis",1,sep="")))
#     RI[iter,2] = rand.index(get(paste("res",2,sep="")),get(paste("resbis",2,sep="")))
#     RI[iter,3] = rand.index(get(paste("res",1,sep="")),get(paste("resbis",2,sep="")))
#     iter = iter + 1
#   }
# }
# 
# fit1 = cmdscale(jaccard_abundance, eig=TRUE, k=k1)
# fit2 = cmdscale(ab_sorensen_abundance, eig=TRUE, k=k1)
# res1 = kmeans(fit1$points,l,nstart = 1000)$cluster
# res2 = kmeans(fit2$points,l,nstart = 1000)$cluster
# rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))

l = 5

# fitJ1 = cmdscale(jaccard_abundance, eig=TRUE, k=2)
# fitJ2 = cmdscale(jaccard_abundance, eig=TRUE, k=3)
# fitJ3 = cmdscale(jaccard_abundance, eig=TRUE, k=4)
# 
# fitS1 = cmdscale(ab_sorensen_abundance, eig=TRUE, k=2)
# fitS2 = cmdscale(ab_sorensen_abundance, eig=TRUE, k=3)
# fitS3 = cmdscale(ab_sorensen_abundance, eig=TRUE, k=4)

fitJ1 = cmdscale(jaccard_abundance, eig=TRUE, k=2)
fitJ2 = cmdscale(jaccard_abundance, eig=TRUE, k=5)
fitJ3 = cmdscale(jaccard_abundance, eig=TRUE, k=10)

fitS1 = cmdscale(ab_sorensen_abundance, eig=TRUE, k=2)
fitS2 = cmdscale(ab_sorensen_abundance, eig=TRUE, k=5)
fitS3 = cmdscale(ab_sorensen_abundance, eig=TRUE, k=10)

fitB1 = cmdscale(braycurtis_abundance, eig=TRUE, k=2)
fitB2 = cmdscale(braycurtis_abundance, eig=TRUE, k=5)
fitB3 = cmdscale(braycurtis_abundance, eig=TRUE, k=10)

resJ1 = kmeans(fitJ1$points,l,nstart = 1000)$cluster
resJ2 = kmeans(fitJ2$points,l,nstart = 1000)$cluster
resJ3 = kmeans(fitJ3$points,l,nstart = 1000)$cluster

resS1 = kmeans(fitS1$points,l,nstart = 1000)$cluster
resS2 = kmeans(fitS2$points,l,nstart = 1000)$cluster
resS3 = kmeans(fitS3$points,l,nstart = 1000)$cluster

resB1 = kmeans(fitB1$points,l,nstart = 1000)$cluster
resB2 = kmeans(fitB2$points,l,nstart = 1000)$cluster
resB3 = kmeans(fitB3$points,l,nstart = 1000)$cluster

RI = matrix(NA,6,6)
RI[1,] = c(rand.index(get(paste("resJ",1,sep="")),get(paste("resJ",1,sep=""))), 
           rand.index(get(paste("resS",1,sep="")),get(paste("resS",1,sep=""))), 
           rand.index(get(paste("resB",1,sep="")),get(paste("resB",1,sep=""))),
           rand.index(get(paste("resJ",1,sep="")),get(paste("resS",1,sep=""))),
           rand.index(get(paste("resJ",1,sep="")),get(paste("resB",1,sep=""))),
           rand.index(get(paste("resS",1,sep="")),get(paste("resB",1,sep=""))))

RI[2,] = c(rand.index(get(paste("resJ",1,sep="")),get(paste("resJ",2,sep=""))), 
           rand.index(get(paste("resS",1,sep="")),get(paste("resS",2,sep=""))), 
           rand.index(get(paste("resB",1,sep="")),get(paste("resB",2,sep=""))),
           rand.index(get(paste("resJ",1,sep="")),get(paste("resS",2,sep=""))),
           rand.index(get(paste("resJ",1,sep="")),get(paste("resB",2,sep=""))),
           rand.index(get(paste("resS",1,sep="")),get(paste("resB",2,sep=""))))

RI[3,] = c(rand.index(get(paste("resJ",1,sep="")),get(paste("resJ",3,sep=""))), 
           rand.index(get(paste("resS",1,sep="")),get(paste("resS",3,sep=""))), 
           rand.index(get(paste("resB",1,sep="")),get(paste("resB",3,sep=""))),
           rand.index(get(paste("resJ",1,sep="")),get(paste("resS",3,sep=""))),
           rand.index(get(paste("resJ",1,sep="")),get(paste("resB",3,sep=""))),
           rand.index(get(paste("resS",1,sep="")),get(paste("resB",3,sep=""))))

RI[4,] = c(rand.index(get(paste("resJ",2,sep="")),get(paste("resJ",2,sep=""))), 
           rand.index(get(paste("resS",2,sep="")),get(paste("resS",2,sep=""))), 
           rand.index(get(paste("resB",2,sep="")),get(paste("resB",2,sep=""))),
           rand.index(get(paste("resJ",2,sep="")),get(paste("resS",2,sep=""))),
           rand.index(get(paste("resJ",2,sep="")),get(paste("resB",2,sep=""))),
           rand.index(get(paste("resS",2,sep="")),get(paste("resB",2,sep=""))))

RI[5,] = c(rand.index(get(paste("resJ",2,sep="")),get(paste("resJ",3,sep=""))), 
           rand.index(get(paste("resS",2,sep="")),get(paste("resS",3,sep=""))), 
           rand.index(get(paste("resB",2,sep="")),get(paste("resB",3,sep=""))),
           rand.index(get(paste("resJ",2,sep="")),get(paste("resS",3,sep=""))),
           rand.index(get(paste("resJ",2,sep="")),get(paste("resB",3,sep=""))),
           rand.index(get(paste("resS",2,sep="")),get(paste("resB",3,sep=""))))

RI[6,] = c(rand.index(get(paste("resJ",3,sep="")),get(paste("resJ",3,sep=""))), 
           rand.index(get(paste("resS",3,sep="")),get(paste("resS",3,sep=""))), 
           rand.index(get(paste("resB",3,sep="")),get(paste("resB",3,sep=""))),
           rand.index(get(paste("resJ",3,sep="")),get(paste("resS",3,sep=""))),
           rand.index(get(paste("resJ",3,sep="")),get(paste("resB",3,sep=""))),
           rand.index(get(paste("resS",3,sep="")),get(paste("resB",3,sep=""))))

###### Prevalence matrices


# ARI1 = matrix(NA,100,50)
# RI1 = matrix(NA,100,50)
# 
# for(k in 1:50)
# {
#   fit1 = cmdscale(kulczynski_prevalence,eig=TRUE, k=k)
#   fit2 = cmdscale(ochiai_prevalence,eig=TRUE, k=k)
#   
#   for(nb_run in 1:100)
#   {
#     res1 = kmeans(fit1$points,l,nstart = 25)$cluster
#     res2 = kmeans(fit2$points,l,nstart = 25)$cluster
#     
#     ARI1[nb_run,k] = adjustedRand(get(paste("res",1,sep="")),get(paste("res",2,sep="")),randMethod="Rand")
#     RI1[nb_run,k] = rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
#   }
# }
# 
# ARI = c()
# RI = c()
# 
# for(k in 1:50){
#   ARI[k] = max(ARI1[,k]) 
#   RI[k] = max(RI1[,k]) 
# }


fitK1 = cmdscale(kulczynski_prevalence, eig=TRUE, k=2)
fitK2 = cmdscale(kulczynski_prevalence, eig=TRUE, k=5)
fitK3 = cmdscale(kulczynski_prevalence, eig=TRUE, k=10)

fitO1 = cmdscale(ochiai_prevalence, eig=TRUE, k=2)
fitO2 = cmdscale(ochiai_prevalence, eig=TRUE, k=5)
fitO3 = cmdscale(ochiai_prevalence, eig=TRUE, k=10)

fitW1 = cmdscale(whittaker_prevalence, eig=TRUE, k=2)
fitW2 = cmdscale(whittaker_prevalence, eig=TRUE, k=5)
fitW3 = cmdscale(whittaker_prevalence, eig=TRUE, k=10)

resK1 = kmeans(fitK1$points,l,nstart = 1000)$cluster
resK2 = kmeans(fitK2$points,l,nstart = 1000)$cluster
resK3 = kmeans(fitK3$points,l,nstart = 1000)$cluster

resO1 = kmeans(fitO1$points,l,nstart = 1000)$cluster
resO2 = kmeans(fitO2$points,l,nstart = 1000)$cluster
resO3 = kmeans(fitO3$points,l,nstart = 1000)$cluster

resW1 = kmeans(fitW1$points,l,nstart = 1000)$cluster
resW2 = kmeans(fitW2$points,l,nstart = 1000)$cluster
resW3 = kmeans(fitW3$points,l,nstart = 1000)$cluster

RI = matrix(NA,6,6)
RI[1,] = c(rand.index(get(paste("resK",1,sep="")),get(paste("resK",1,sep=""))), 
           rand.index(get(paste("resO",1,sep="")),get(paste("resO",1,sep=""))), 
           rand.index(get(paste("resW",1,sep="")),get(paste("resW",1,sep=""))),
           rand.index(get(paste("resK",1,sep="")),get(paste("resO",1,sep=""))),
           rand.index(get(paste("resK",1,sep="")),get(paste("resW",1,sep=""))),
           rand.index(get(paste("resO",1,sep="")),get(paste("resW",1,sep=""))))

RI[2,] = c(rand.index(get(paste("resK",1,sep="")),get(paste("resK",2,sep=""))), 
           rand.index(get(paste("resO",1,sep="")),get(paste("resO",2,sep=""))), 
           rand.index(get(paste("resW",1,sep="")),get(paste("resW",2,sep=""))),
           rand.index(get(paste("resK",1,sep="")),get(paste("resO",2,sep=""))),
           rand.index(get(paste("resK",1,sep="")),get(paste("resW",2,sep=""))),
           rand.index(get(paste("resO",1,sep="")),get(paste("resW",2,sep=""))))

RI[3,] = c(rand.index(get(paste("resK",1,sep="")),get(paste("resK",3,sep=""))), 
           rand.index(get(paste("resO",1,sep="")),get(paste("resO",3,sep=""))), 
           rand.index(get(paste("resW",1,sep="")),get(paste("resW",3,sep=""))),
           rand.index(get(paste("resK",1,sep="")),get(paste("resO",3,sep=""))),
           rand.index(get(paste("resK",1,sep="")),get(paste("resW",3,sep=""))),
           rand.index(get(paste("resO",1,sep="")),get(paste("resW",3,sep=""))))

RI[4,] = c(rand.index(get(paste("resK",2,sep="")),get(paste("resK",2,sep=""))), 
           rand.index(get(paste("resO",2,sep="")),get(paste("resO",2,sep=""))), 
           rand.index(get(paste("resW",2,sep="")),get(paste("resW",2,sep=""))),
           rand.index(get(paste("resK",2,sep="")),get(paste("resO",2,sep=""))),
           rand.index(get(paste("resK",2,sep="")),get(paste("resW",2,sep=""))),
           rand.index(get(paste("resO",2,sep="")),get(paste("resW",2,sep=""))))

RI[5,] = c(rand.index(get(paste("resK",2,sep="")),get(paste("resK",3,sep=""))), 
           rand.index(get(paste("resO",2,sep="")),get(paste("resO",3,sep=""))), 
           rand.index(get(paste("resW",2,sep="")),get(paste("resW",3,sep=""))),
           rand.index(get(paste("resK",2,sep="")),get(paste("resO",3,sep=""))),
           rand.index(get(paste("resK",2,sep="")),get(paste("resW",3,sep=""))),
           rand.index(get(paste("resO",2,sep="")),get(paste("resW",3,sep=""))))

RI[6,] = c(rand.index(get(paste("resK",3,sep="")),get(paste("resK",3,sep=""))), 
           rand.index(get(paste("resO",3,sep="")),get(paste("resO",3,sep=""))), 
           rand.index(get(paste("resW",3,sep="")),get(paste("resW",3,sep=""))),
           rand.index(get(paste("resK",3,sep="")),get(paste("resO",3,sep=""))),
           rand.index(get(paste("resK",3,sep="")),get(paste("resW",3,sep=""))),
           rand.index(get(paste("resO",3,sep="")),get(paste("resW",3,sep=""))))


############################## k-NNG + spectral clustering


###### Abundance matrices

# l = 5
# 
# ARI1 = matrix(NA,100,50)
# RI1 = matrix(NA,100,50)
# 
# for(k in 1:50)
# {
#   kNN1 =  make.kNNG(jaccard_abundance, k = k, symm = TRUE, weight = FALSE)
#   kNN2 = make.kNNG(ab_sorensen_abundance, k = k, symm = TRUE, weight = FALSE)
#   
#   for(nb_run in 1:100)
#   {
#     res1 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#     res2 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#     
#     ARI1[nb_run,k] = adjustedRand(get(paste("res",1,sep="")),get(paste("res",2,sep="")),randMethod="Rand")
#     RI1[nb_run,k] = rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
#   }
# }
# 
# ARI = c()
# RI = c()
# 
# for(k in 1:50){
#   ARI[k] = max(ARI1[,k]) 
#   RI[k] = max(RI1[,k]) 
# }

fitJ1 = make.kNNG(jaccard_abundance, k = 10, symm = TRUE, weight = FALSE)
fitJ2 = make.kNNG(jaccard_abundance, k = 25, symm = TRUE, weight = FALSE)
fitJ3 = make.kNNG(jaccard_abundance, k = 40, symm = TRUE, weight = FALSE)

fitS1 = make.kNNG(ab_sorensen_abundance, k = 10, symm = TRUE, weight = FALSE)
fitS2 = make.kNNG(ab_sorensen_abundance, k = 25, symm = TRUE, weight = FALSE)
fitS3 = make.kNNG(ab_sorensen_abundance, k = 40, symm = TRUE, weight = FALSE)

fitB1 = make.kNNG(braycurtis_abundance, k = 10, symm = TRUE, weight = FALSE)
fitB2 = make.kNNG(braycurtis_abundance, k = 25, symm = TRUE, weight = FALSE)
fitB3 = make.kNNG(braycurtis_abundance, k = 40, symm = TRUE, weight = FALSE)

RI = matrix(0,6,6)

for(iter in 1:100){
  resJ1 = spectral.clustering(fitJ1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  resJ2 = spectral.clustering(fitJ2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  resJ3 = spectral.clustering(fitJ3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  resS1 = spectral.clustering(fitS1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  resS2 = spectral.clustering(fitS2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  resS3 = spectral.clustering(fitS3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  resB1 = spectral.clustering(fitB1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  resB2 = spectral.clustering(fitB2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  resB3 = spectral.clustering(fitB3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  RI[1,] = RI[1,] + c(rand.index(get(paste("resJ",1,sep="")),get(paste("resJ",1,sep=""))), 
                      rand.index(get(paste("resS",1,sep="")),get(paste("resS",1,sep=""))), 
                      rand.index(get(paste("resB",1,sep="")),get(paste("resB",1,sep=""))),
                      rand.index(get(paste("resJ",1,sep="")),get(paste("resS",1,sep=""))),
                      rand.index(get(paste("resJ",1,sep="")),get(paste("resB",1,sep=""))),
                      rand.index(get(paste("resS",1,sep="")),get(paste("resB",1,sep=""))))
  
  RI[2,] = RI[2,] + c(rand.index(get(paste("resJ",1,sep="")),get(paste("resJ",2,sep=""))), 
                      rand.index(get(paste("resS",1,sep="")),get(paste("resS",2,sep=""))), 
                      rand.index(get(paste("resB",1,sep="")),get(paste("resB",2,sep=""))),
                      rand.index(get(paste("resJ",1,sep="")),get(paste("resS",2,sep=""))),
                      rand.index(get(paste("resJ",1,sep="")),get(paste("resB",2,sep=""))),
                      rand.index(get(paste("resS",1,sep="")),get(paste("resB",2,sep=""))))
  
  RI[3,] = RI[3,] + c(rand.index(get(paste("resJ",1,sep="")),get(paste("resJ",3,sep=""))), 
                      rand.index(get(paste("resS",1,sep="")),get(paste("resS",3,sep=""))), 
                      rand.index(get(paste("resB",1,sep="")),get(paste("resB",3,sep=""))),
                      rand.index(get(paste("resJ",1,sep="")),get(paste("resS",3,sep=""))),
                      rand.index(get(paste("resJ",1,sep="")),get(paste("resB",3,sep=""))),
                      rand.index(get(paste("resS",1,sep="")),get(paste("resB",3,sep=""))))
  
  RI[4,] = RI[4,] + c(rand.index(get(paste("resJ",2,sep="")),get(paste("resJ",2,sep=""))), 
                      rand.index(get(paste("resS",2,sep="")),get(paste("resS",2,sep=""))), 
                      rand.index(get(paste("resB",2,sep="")),get(paste("resB",2,sep=""))),
                      rand.index(get(paste("resJ",2,sep="")),get(paste("resS",2,sep=""))),
                      rand.index(get(paste("resJ",2,sep="")),get(paste("resB",2,sep=""))),
                      rand.index(get(paste("resS",2,sep="")),get(paste("resB",2,sep=""))))
  
  RI[5,] = RI[5,] + c(rand.index(get(paste("resJ",2,sep="")),get(paste("resJ",3,sep=""))), 
                      rand.index(get(paste("resS",2,sep="")),get(paste("resS",3,sep=""))), 
                      rand.index(get(paste("resB",2,sep="")),get(paste("resB",3,sep=""))),
                      rand.index(get(paste("resJ",2,sep="")),get(paste("resS",3,sep=""))),
                      rand.index(get(paste("resJ",2,sep="")),get(paste("resB",3,sep=""))),
                      rand.index(get(paste("resS",2,sep="")),get(paste("resB",3,sep=""))))
  
  RI[6,] = RI[6,] + c(rand.index(get(paste("resJ",3,sep="")),get(paste("resJ",3,sep=""))), 
                      rand.index(get(paste("resS",3,sep="")),get(paste("resS",3,sep=""))), 
                      rand.index(get(paste("resB",3,sep="")),get(paste("resB",3,sep=""))),
                      rand.index(get(paste("resJ",3,sep="")),get(paste("resS",3,sep=""))),
                      rand.index(get(paste("resJ",3,sep="")),get(paste("resB",3,sep=""))),
                      rand.index(get(paste("resS",3,sep="")),get(paste("resB",3,sep=""))))
}

RI = RI/100

###### Prevalence matrices


# ARI1 = matrix(NA,100,50)
# RI1 = matrix(NA,100,50)
#
# for(k in 1:50)
# {
#   kNN1 =  make.kNNG(ochiai_prevalence, k = k, symm = TRUE, weight = FALSE)
#   kNN2 = make.kNNG(kulczynski_prevalence, k = k, symm = TRUE, weight = FALSE)
#   
#   for(nb_run in 1:100)
#   {
#     res1 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#     res2 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#     
#     ARI1[nb_run,k] = adjustedRand(get(paste("res",1,sep="")),get(paste("res",2,sep="")),randMethod="Rand")
#     RI1[nb_run,k] = rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
#   }
# }
# 
# ARI = c()
# RI = c()
# 
# for(k in 1:50){
#   ARI[k] = max(ARI1[,k]) 
#   RI[k] = max(RI1[,k]) 
# }

fitK1 = make.kNNG(kulczynski_prevalence, k = 10, symm = TRUE, weight = FALSE)
fitK2 = make.kNNG(kulczynski_prevalence, k = 25, symm = TRUE, weight = FALSE)
fitK3 = make.kNNG(kulczynski_prevalence, k = 40, symm = TRUE, weight = FALSE)

fitO1 = make.kNNG(ochiai_prevalence, k = 10, symm = TRUE, weight = FALSE)
fitO2 = make.kNNG(ochiai_prevalence, k = 25, symm = TRUE, weight = FALSE)
fitO3 = make.kNNG(ochiai_prevalence, k = 40, symm = TRUE, weight = FALSE)

fitW1 = make.kNNG(whittaker_prevalence, k = 10, symm = TRUE, weight = FALSE)
fitW2 = make.kNNG(whittaker_prevalence, k = 25, symm = TRUE, weight = FALSE)
fitW3 = make.kNNG(whittaker_prevalence, k = 40, symm = TRUE, weight = FALSE)

RI = matrix(0,6,6)

for(iter in 1:100){
  resK1 = spectral.clustering(fitK1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  resK2 = spectral.clustering(fitK2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  resK3 = spectral.clustering(fitK3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  resO1 = spectral.clustering(fitO1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  resO2 = spectral.clustering(fitO2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  resO3 = spectral.clustering(fitO3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  resW1 = spectral.clustering(fitW1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  resW2 = spectral.clustering(fitW2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  resW3 = spectral.clustering(fitW3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  RI[1,] = RI[1,] + c(rand.index(get(paste("resK",1,sep="")),get(paste("resK",1,sep=""))), 
                      rand.index(get(paste("resO",1,sep="")),get(paste("resO",1,sep=""))), 
                      rand.index(get(paste("resW",1,sep="")),get(paste("resW",1,sep=""))),
                      rand.index(get(paste("resK",1,sep="")),get(paste("resO",1,sep=""))),
                      rand.index(get(paste("resK",1,sep="")),get(paste("resW",1,sep=""))),
                      rand.index(get(paste("resO",1,sep="")),get(paste("resW",1,sep=""))))
  
  RI[2,] = RI[2,] + c(rand.index(get(paste("resK",1,sep="")),get(paste("resK",2,sep=""))), 
                      rand.index(get(paste("resO",1,sep="")),get(paste("resO",2,sep=""))), 
                      rand.index(get(paste("resW",1,sep="")),get(paste("resW",2,sep=""))),
                      rand.index(get(paste("resK",1,sep="")),get(paste("resO",2,sep=""))),
                      rand.index(get(paste("resK",1,sep="")),get(paste("resW",2,sep=""))),
                      rand.index(get(paste("resO",1,sep="")),get(paste("resW",2,sep=""))))
  
  RI[3,] = RI[3,] + c(rand.index(get(paste("resK",1,sep="")),get(paste("resK",3,sep=""))), 
                      rand.index(get(paste("resO",1,sep="")),get(paste("resO",3,sep=""))), 
                      rand.index(get(paste("resW",1,sep="")),get(paste("resW",3,sep=""))),
                      rand.index(get(paste("resK",1,sep="")),get(paste("resO",3,sep=""))),
                      rand.index(get(paste("resK",1,sep="")),get(paste("resW",3,sep=""))),
                      rand.index(get(paste("resO",1,sep="")),get(paste("resW",3,sep=""))))
  
  RI[4,] = RI[4,] + c(rand.index(get(paste("resK",2,sep="")),get(paste("resK",2,sep=""))), 
                      rand.index(get(paste("resO",2,sep="")),get(paste("resO",2,sep=""))), 
                      rand.index(get(paste("resW",2,sep="")),get(paste("resW",2,sep=""))),
                      rand.index(get(paste("resK",2,sep="")),get(paste("resO",2,sep=""))),
                      rand.index(get(paste("resK",2,sep="")),get(paste("resW",2,sep=""))),
                      rand.index(get(paste("resO",2,sep="")),get(paste("resW",2,sep=""))))
  
  
  RI[5,] = RI[5,] + c(rand.index(get(paste("resK",2,sep="")),get(paste("resK",3,sep=""))), 
                      rand.index(get(paste("resO",2,sep="")),get(paste("resO",3,sep=""))), 
                      rand.index(get(paste("resW",2,sep="")),get(paste("resW",3,sep=""))),
                      rand.index(get(paste("resK",2,sep="")),get(paste("resO",3,sep=""))),
                      rand.index(get(paste("resK",2,sep="")),get(paste("resW",3,sep=""))),
                      rand.index(get(paste("resO",2,sep="")),get(paste("resW",3,sep=""))))
  
  RI[6,] = RI[6,] + c(rand.index(get(paste("resK",3,sep="")),get(paste("resK",3,sep=""))), 
                      rand.index(get(paste("resO",3,sep="")),get(paste("resO",3,sep=""))), 
                      rand.index(get(paste("resW",3,sep="")),get(paste("resW",3,sep=""))),
                      rand.index(get(paste("resK",3,sep="")),get(paste("resO",3,sep=""))),
                      rand.index(get(paste("resK",3,sep="")),get(paste("resW",3,sep=""))),
                      rand.index(get(paste("resO",3,sep="")),get(paste("resW",3,sep=""))))
  
}

RI = RI/100


##################################### Comparison between clusterings 

######### Abundance matrices

fit1 = cmdscale(jaccard_abundance, eig=TRUE, k=10)
fit2 = cmdscale(ab_sorensen_abundance,eig=TRUE, k=10)
res1 = kmeans(fit1$points,l)$cluster
res2 = kmeans(fit2$points,l)$cluster

kNN1 =  make.kNNG(jaccard_abundance, k = 41, symm = TRUE, weight = FALSE)
kNN2 = make.kNNG(ab_sorensen_abundance, k = 41, symm = TRUE, weight = FALSE)
res3 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res4 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
rand.index(get(paste("res",3,sep="")),get(paste("res",4,sep="")))
rand.index(get(paste("res",1,sep="")),get(paste("res",3,sep="")))
rand.index(get(paste("res",2,sep="")),get(paste("res",4,sep="")))

######### Prevalence matrices

fit1 = cmdscale(ochiai_prevalence, eig=TRUE, k=10)
fit2 = cmdscale(kulczynski_prevalence,eig=TRUE, k=10)
res1 = kmeans(fit1$points,l)$cluster
res2 = kmeans(fit2$points,l)$cluster

kNN1 =  make.kNNG(ochiai_prevalence, k = 35, symm = TRUE, weight = FALSE)
kNN2 = make.kNNG(kulczynski_prevalence, k = 35, symm = TRUE, weight = FALSE)
res3 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res4 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
rand.index(get(paste("res",3,sep="")),get(paste("res",4,sep="")))
rand.index(get(paste("res",1,sep="")),get(paste("res",3,sep="")))
rand.index(get(paste("res",2,sep="")),get(paste("res",4,sep="")))

# fit1 = cmdscale(ochiai_prevalence, eig=TRUE, k=10)
# fit2 = cmdscale(ochiai_prevalence,eig=TRUE, k=4)
# fit3 = cmdscale(kulczynski_prevalence, eig=TRUE, k=10)
# fit4 = cmdscale(kulczynski_prevalence,eig=TRUE, k=4)
# res1 = kmeans(fit1$points,l)$cluster
# res2 = kmeans(fit2$points,l)$cluster
# res3 = kmeans(fit3$points,l)$cluster
# res4 = kmeans(fit4$points,l)$cluster
# rand.index(get(paste("res",1,sep="")),get(paste("res",3,sep="")))
# rand.index(get(paste("res",2,sep="")),get(paste("res",4,sep="")))
# rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
# rand.index(get(paste("res",3,sep="")),get(paste("res",4,sep="")))


############################### To search which k maximises similarity between graphs

############# Abundance matrices

res = c()

for(k in 1:50){
  
  kNN1 = make.kNNG(jaccard_abundance, k = k, symm = TRUE, weight = FALSE)

  kNN2 = make.kNNG(ab_jaccard_abundance, k = k, symm = TRUE, weight = FALSE)

  kNN3 = make.kNNG(braycurtis_abundance, k = k, symm = TRUE, weight = FALSE)

  kNN4 = make.kNNG(ab_ochiai_abundance, k = k, symm = TRUE, weight = FALSE)

  kNN5 = make.kNNG(ab_sorensen_abundance, k = k, symm = TRUE, weight = FALSE)

  kNN6 = make.kNNG(simka_jaccard_abundance, k = k, symm = TRUE, weight = FALSE)
  
  index = 0
  
  for(j in 1:5){
    for(m in (j+1):6){
      index = index + GARI(get(paste("kNN",j,sep="")),get(paste("kNN",m,sep="")))
    }
  }
  res[k] = index/15
}

res
which.max(res)
max(res)

############# Prevalence matrices

res = c()

for(k in 1:50){
  
  kNN7 = make.kNNG(chord_prevalence, k = k, symm = TRUE, weight = FALSE)
  
  kNN8 = make.kNNG(jaccard_prevalence, k = k, symm = TRUE, weight = FALSE)
  
  kNN9 = make.kNNG(kulczynski_prevalence, k = k, symm = TRUE, weight = FALSE)
  
  kNN10 = make.kNNG(ochiai_prevalence, k = k, symm = TRUE, weight = FALSE)
  
  kNN11 = make.kNNG(whittaker_prevalence, k = k, symm = TRUE, weight = FALSE)
  
  kNN12 = make.kNNG(simka_jaccard_prevalence, k = k, symm = TRUE, weight = FALSE)
  
  index = 0
  
  for(j in 7:11){
    for(m in (j+1):12){
      index = index + GARI(get(paste("kNN",j,sep="")),get(paste("kNN",m,sep="")))
    }
  }
  res[k] = index/15
}

res
which.max(res)
max(res)

####################################################################### Multivariate analysis 

rm(list=objects())

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)


############################################################ Source custom scripts

source('~/metaord/utils.R')

############################################################ Design

data.wd <- "/home/rnarci/Bureau/CDDrnarci/Donnees/"
design = read.table(file=file.path(data.wd, "param_bioadvection.csv"),sep="",header=TRUE)
rownames(design) <- design$Sample

############################################################ Import longitude and latitude of the stations

gps = read.table(file=file.path(data.wd,"GPScoordinates2.csv"),header=TRUE,sep="")
rownames(gps) <- gps$Station

############################################################ Import data

size_fraction1 = "0-0.2"
size_fraction2 = "0.22-3"
size_fraction3 = "0.8-5"
size_fraction4 = "5-20"
size_fraction5 = "20-180"
size_fraction6 = "180-2000"

import_data(size_fraction3, samples = rownames(design))

############################################################ Subset design and (longitude,latitude)

design <- design[metagenomic_sample, ]
design <- design[,-c(1,6,7)]
# design <- design[,-c(1,5,6,7,13)]

gps <- gps[metagenomic_sample,]
gps <- gps[,-1]

library(mice)
mice = complete(mice(design,method="norm.predict",m=1))
design$Phosphates = mice$Phosphates
design$NO2NO3 = mice$NO2NO3

lagrangien = read.table(file=file.path(data.wd,"tarrive_min_surface_1000.csv"),header=TRUE, sep="")
colnames(lagrangien) <- rownames(lagrangien)
lagrangien = lagrangien[metagenomic_sample,metagenomic_sample]

######################################## To calculate predicted distances

D = jaccard_abundance
n = dim(D)[1]

parms1 = c(design$Temperature,design$SSD,design$SI_temperature,design$Depth,design$Silicate,design$SI_nitrates)
X1 = matrix(parms1,nrow=n,ncol=6)
parms2 = design$Temperature
X2 = matrix(parms2,nrow=n,ncol=1)

H1 = X1%*%solve(t(X1)%*%X1)%*%t(X1)
H2 = X2%*%solve(t(X2)%*%X2)%*%t(X2)

D_without_MDS1 = matrix(NA,n,n)
D_without_MDS2 = matrix(NA,n,n)

J = diag(rep(1,n)) - matrix(1,n,n)/n
G = -0.5*J%*%D^2%*%J

for(i in 1:n){
  for(j in 1:n){
    D_without_MDS1[i,j] = sqrt(H1[i,]%*%G%*%H1[,i] + H1[j,]%*%G%*%H1[,j] - 2*H1[i,]%*%G%*%H1[,j])
    D_without_MDS2[i,j] = sqrt(H2[i,]%*%G%*%H2[,i] + H2[j,]%*%G%*%H2[,j] - 2*H2[i,]%*%G%*%H2[,j])
  }
}

fit1 = cmdscale(D, eig=TRUE, k=10)
fit2 = cmdscale(D_without_MDS1, eig=TRUE, k=10)
fit3 = cmdscale(D_without_MDS2, eig=TRUE, k=10)

library(loe)
library(fcd)
library(clues)

kNN1 = make.kNNG(D, k = 37, symm = TRUE, weight = FALSE)
kNN2 = make.kNNG(D_without_MDS1, k = 37, symm = TRUE, weight = FALSE)
kNN3 = make.kNNG(D_without_MDS2, k = 37, symm = TRUE, weight = FALSE)

l = 3

res1 = kmeans(fit1$points,l, nstart = 1000)$cluster
res2 = kmeans(fit2$points,l, nstart = 1000)$cluster
res3 = kmeans(fit3$points,l, nstart = 1000)$cluster

res4 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res5 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res6 = spectral.clustering(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)


########## Boxplots : clustering on numeric distances

temperature_data = data.frame(clus.raw = res1,
                              clus.predicted.all.variables = res2
                              ) %>% 
  mutate(Temperature = design$Temperature)

tdata <- gather(temperature_data, key = "Clustering", value = "Group", clus.raw:clus.predicted.all.variables)

ggplot(tdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 2)

########## Boxplots : clustering on ordinal distances

# temperature_data = data.frame(clus.raw = res4,
#                               clus.predicted.all.variables = res5,
#                               clus.predicted.only.temperature = res6) %>% 
#   mutate(Temperature = design$Temperature)
# 
# tdata <- gather(temperature_data, key = "Clustering", value = "Group", clus.raw:clus.predicted.only.temperature)
# 
# ggplot(tdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 3)

################################## Map resulting clusters on raw distances (ordinal setting)

library(rworldmap)
newmap <- getMap(resolution = "li")

color = c()
labels = c()
cl = res1

for(i in 1:length(cl)){
  if(cl[i]==1){
    color[i] = "coral4"
  }
  if(cl[i]==2){
    color[i] = "darkviolet"
  }
  if(cl[i]==3){
    color[i] = "black"
  }
  if(str_length(rownames(gps)[i]) == 5){
    labels[i] = substr(rownames(gps)[i],1,1)
  }
  if(str_length(rownames(gps)[i]) == 6){
    labels[i] = substr(rownames(gps)[i],1,2)
  }
  if(str_length(rownames(gps)[i]) == 7){
    labels[i] = substr(rownames(gps)[i],1,3)
  }
}

DCM_indices = which(str_detect(rownames(gps),"DCM")==TRUE)
SUR_indices = which(str_detect(rownames(gps),"SUR")==TRUE)

SUR_only_indices = c()
a = 1
for(i in 1:length(cl)){
  if(length(which(labels[SUR_indices][i]!=labels[DCM_indices]))==length(labels[DCM_indices])){
    SUR_only_indices[a] = SUR_indices[i]  
    a = a + 1
  }
}

plot(newmap, xlim = c(-180, 90), ylim = c(-75, 75), asp = 1)
points(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], col = color[DCM_indices], cex = 1, pch = 17)
points(gps$Mean_longitude[SUR_indices], gps$Mean_latitude[SUR_indices] + 2, col = color[SUR_indices], cex = 1, pch = 16)
text(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], labels = labels[DCM_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)
text(gps$Mean_longitude[SUR_only_indices], gps$Mean_latitude[SUR_only_indices] + 2, labels = labels[SUR_only_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)



########## Map of resulting clusters : raw vs predict. with Temperature only

color = c()
labels = c()
cl = res3

for(i in 1:length(cl)){
  if(cl[i]==1){
    color[i] = "blue"
  }
  if(cl[i]==2){
    color[i] = "red"
  }
  if(cl[i]==3){
    color[i] = "gold4"
  }
  if(str_length(rownames(gps)[i]) == 5){
    labels[i] = substr(rownames(gps)[i],1,1)
  }
  if(str_length(rownames(gps)[i]) == 6){
    labels[i] = substr(rownames(gps)[i],1,2)
  }
  if(str_length(rownames(gps)[i]) == 7){
    labels[i] = substr(rownames(gps)[i],1,3)
  }
}

DCM_indices = which(str_detect(rownames(gps),"DCM")==TRUE)
SUR_indices = which(str_detect(rownames(gps),"SUR")==TRUE)

SUR_only_indices = c()
a = 1
for(i in 1:length(cl)){
  if(length(which(labels[SUR_indices][i]!=labels[DCM_indices]))==length(labels[DCM_indices])){
    SUR_only_indices[a] = SUR_indices[i]  
    a = a + 1
  }
}

plot(newmap, xlim = c(-180, 90), ylim = c(-75, 75), asp = 1)
points(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], col = color[DCM_indices], cex = 1, pch = 17)
points(gps$Mean_longitude[SUR_indices], gps$Mean_latitude[SUR_indices] + 2, col = color[SUR_indices], cex = 1, pch = 16)
text(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], labels = labels[DCM_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)
text(gps$Mean_longitude[SUR_only_indices], gps$Mean_latitude[SUR_only_indices] + 2, labels = labels[SUR_only_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)




########## Map of resulting clusters : raw vs predict. with all covariates

color = c()
labels = c()
cl = res2

for(i in 1:length(cl)){
  if(cl[i]==1){
    color[i] = "blue"
  }
  if(cl[i]==2){
    color[i] = "red"
  }
  if(cl[i]==3){
    color[i] = "gold4"
  }
  if(str_length(rownames(gps)[i]) == 5){
    labels[i] = substr(rownames(gps)[i],1,1)
  }
  if(str_length(rownames(gps)[i]) == 6){
    labels[i] = substr(rownames(gps)[i],1,2)
  }
  if(str_length(rownames(gps)[i]) == 7){
    labels[i] = substr(rownames(gps)[i],1,3)
  }
}

DCM_indices = which(str_detect(rownames(gps),"DCM")==TRUE)
SUR_indices = which(str_detect(rownames(gps),"SUR")==TRUE)

SUR_only_indices = c()
a = 1
for(i in 1:length(cl)){
  if(length(which(labels[SUR_indices][i]!=labels[DCM_indices]))==length(labels[DCM_indices])){
    SUR_only_indices[a] = SUR_indices[i]  
    a = a + 1
  }
}

plot(newmap, xlim = c(-180, 90), ylim = c(-75, 75), asp = 1)
points(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], col = color[DCM_indices], cex = 1, pch = 17)
points(gps$Mean_longitude[SUR_indices], gps$Mean_latitude[SUR_indices] + 2, col = color[SUR_indices], cex = 1, pch = 16)
text(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], labels = labels[DCM_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)
text(gps$Mean_longitude[SUR_only_indices], gps$Mean_latitude[SUR_only_indices] + 2, labels = labels[SUR_only_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)

########## Plot a distance matrix as a heatmap

library(ggplot2)
library(reshape2)

plot_dist_as_heatmap <- function(dist, order = NULL, title = NULL,
                                 low = "#B1F756", high = "#132B13") {
  data <- melt(as(dist, "matrix"))
  colnames(data) <- c("x", "y", "distance")
  if (!is.null(order)) {
    data$x <- factor(data$x, levels = order)
    data$y <- factor(data$y, levels = order)
  }
  p <- ggplot(data, aes(x = x, y = y, fill = distance)) + geom_tile()
  p <- p + theme(axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.text.y = element_blank()) +
    scale_fill_gradient(low = low, high = high)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

####### Remove colnames and rownames

D_bis = matrix(NA,n,n)

for(i in 1:n){
  D_bis[i,] = D[i,]
}

plot_dist_as_heatmap(as.dist(D_bis), order = order(as.numeric(labels)), title = "Heatmap for raw distance matrix")
plot_dist_as_heatmap(as.dist(D_without_MDS2), order = order(as.numeric(labels)), title = "Heatmap for predicted distance matrix (only temperature)")
plot_dist_as_heatmap(as.dist(D_without_MDS1), order = order(as.numeric(labels)), title = "Heatmap for predicted distance matrix (all variables)")

######## Plot of two distance matrices

plot(as.dist(D),as.dist(D_without_MDS1),xlab="Distances brutes",ylab="Distances prédites",pch=16,col="blue", main="Toutes les covariables",cex.main=1.8,cex.lab=1.3)
plot(as.dist(D),as.dist(D_without_MDS2),xlab="Distances brutes",ylab="Distances prédites",pch=16,col="blue", main="Temperature uniquement",cex.main=1.8,cex.lab=1.3)
plot(as.dist(D),as.dist(lagrangien),xlab="Distances brutes",ylab="Distances Lagrangiennes",pch=16,col="blue",cex.lab=1.3)

######################################## To calculate corrected distances

D_without_MDS1 = matrix(NA,n,n)
D_without_MDS2 = matrix(NA,n,n)
I = diag(rep(1,n))

for(i in 1:n){
  for(j in 1:n){
    D_without_MDS1[i,j] = sqrt((I-H1)[i,]%*%G%*%(I-H1)[,i] + (I-H1)[j,]%*%G%*%(I-H1)[,j] - 2*(I-H1)[i,]%*%G%*%(I-H1)[,j])
    D_without_MDS2[i,j] = sqrt((I-H2)[i,]%*%G%*%(I-H2)[,i] + (I-H2)[j,]%*%G%*%(I-H2)[,j] - 2*(I-H2)[i,]%*%G%*%(I-H2)[,j])
  }
}

fit1 = cmdscale(D, eig=TRUE, k=10)
fit2 = cmdscale(D_without_MDS1, eig=TRUE, k=10)
fit3 = cmdscale(D_without_MDS2, eig=TRUE, k=10)

kNN1 = make.kNNG(D, k = 37, symm = TRUE, weight = FALSE)
kNN2 = make.kNNG(D_without_MDS1, k = 37, symm = TRUE, weight = FALSE)
kNN3 = make.kNNG(D_without_MDS2, k = 37, symm = TRUE, weight = FALSE)

res1 = kmeans(fit1$points,l, nstart = 1000)$cluster
res2 = kmeans(fit2$points,l, nstart = 1000)$cluster
res3 = kmeans(fit3$points,l, nstart = 1000)$cluster


res4 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res5 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res6 = spectral.clustering(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)


########## Boxplots : clustering on numeric distances

temperature_data = data.frame(clus.raw = res1,
                              clus.corrected.all.variables = res2) %>% 
  mutate(Temperature = design$Temperature)

tdata <- gather(temperature_data, key = "Clustering", value = "Group", clus.raw:clus.corrected.all.variables)

ggplot(tdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 2)


########## Boxplots : clustering on ordinal distances

# temperature_data = data.frame(clus.raw = res4,
#                               clus.corrected.all.variables = res5,
#                               clus.corrected.only.temperature = res6) %>% 
#   mutate(Temperature = design$Temperature)
# 
# tdata <- gather(temperature_data, key = "Clustering", value = "Group", clus.raw:clus.corrected.only.temperature)
# 
# ggplot(tdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 3)


########## Map of resulting clusters : raw vs correct. with Temperature only

color = c()
labels = c()
cl = res3

for(i in 1:length(cl)){
  if(cl[i]==1){
    color[i] = "blue"
  }
  if(cl[i]==2){
    color[i] = "red"
  }
  if(cl[i]==3){
    color[i] = "gold4"
  }
  if(str_length(rownames(gps)[i]) == 5){
    labels[i] = substr(rownames(gps)[i],1,1)
  }
  if(str_length(rownames(gps)[i]) == 6){
    labels[i] = substr(rownames(gps)[i],1,2)
  }
  if(str_length(rownames(gps)[i]) == 7){
    labels[i] = substr(rownames(gps)[i],1,3)
  }
}

DCM_indices = which(str_detect(rownames(gps),"DCM")==TRUE)
SUR_indices = which(str_detect(rownames(gps),"SUR")==TRUE)

SUR_only_indices = c()
a = 1
for(i in 1:length(cl)){
  if(length(which(labels[SUR_indices][i]!=labels[DCM_indices]))==length(labels[DCM_indices])){
    SUR_only_indices[a] = SUR_indices[i]  
    a = a + 1
  }
}

plot(newmap, xlim = c(-180, 90), ylim = c(-75, 75), asp = 1)
points(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], col = color[DCM_indices], cex = 1, pch = 17)
points(gps$Mean_longitude[SUR_indices], gps$Mean_latitude[SUR_indices] + 2, col = color[SUR_indices], cex = 1, pch = 16)
text(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], labels = labels[DCM_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)
text(gps$Mean_longitude[SUR_only_indices], gps$Mean_latitude[SUR_only_indices] + 2, labels = labels[SUR_only_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)




########## Map of resulting clusters : raw vs correct. with all covariates

color = c()
labels = c()
cl = res2

for(i in 1:length(cl)){
  if(cl[i]==1){
    color[i] = "red"
  }
  if(cl[i]==2){
    color[i] = "gold4"
  }
  if(cl[i]==3){
    color[i] = "blue"
  }
  if(str_length(rownames(gps)[i]) == 5){
    labels[i] = substr(rownames(gps)[i],1,1)
  }
  if(str_length(rownames(gps)[i]) == 6){
    labels[i] = substr(rownames(gps)[i],1,2)
  }
  if(str_length(rownames(gps)[i]) == 7){
    labels[i] = substr(rownames(gps)[i],1,3)
  }
}

DCM_indices = which(str_detect(rownames(gps),"DCM")==TRUE)
SUR_indices = which(str_detect(rownames(gps),"SUR")==TRUE)

SUR_only_indices = c()
a = 1
for(i in 1:length(cl)){
  if(length(which(labels[SUR_indices][i]!=labels[DCM_indices]))==length(labels[DCM_indices])){
    SUR_only_indices[a] = SUR_indices[i]  
    a = a + 1
  }
}

plot(newmap, xlim = c(-180, 90), ylim = c(-75, 75), asp = 1)
points(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], col = color[DCM_indices], cex = 1, pch = 17)
points(gps$Mean_longitude[SUR_indices], gps$Mean_latitude[SUR_indices] + 2, col = color[SUR_indices], cex = 1, pch = 16)
text(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], labels = labels[DCM_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)
text(gps$Mean_longitude[SUR_only_indices], gps$Mean_latitude[SUR_only_indices] + 2, labels = labels[SUR_only_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)

########## Plot a distance matrix as a heatmap

plot_dist_as_heatmap(as.dist(D_without_MDS2), order = order(as.numeric(labels)), title = "Heatmap for corrected distance matrix (only temperature)")
plot_dist_as_heatmap(as.dist(D_without_MDS1), order = order(as.numeric(labels)), title = "Heatmap for corrected distance matrix (all variables)")

######## Plot of two distance matrices

plot(as.dist(D),as.dist(D_without_MDS1),xlab="Distances brutes",ylab="Distances corrigées",pch=16,col="blue",main="Toutes les covariables",cex.main=1.8,cex.lab=1.3)
plot(as.dist(D),as.dist(D_without_MDS2),xlab="Distances brutes",ylab="Distances corrigées",pch=16,col="blue",main="Temperature uniquement",cex.main=1.8,cex.lab=1.3)
