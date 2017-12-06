rm(list=objects())

library(fcd)
library(loe)
library(clues)
library(cccd)
library(fossil)

jaccard_abundance = read.table(file="mat_abundance_jaccard.csv",sep=",")
ab_jaccard_abundance = read.table(file="mat_abundance_ab-jaccard.csv",sep=",")
braycurtis_abundance = read.table(file="mat_abundance_braycurtis.csv",sep=",")
ab_ochiai_abundance = read.table(file="mat_abundance_ab-ochiai.csv",sep=",")
ab_sorensen_abundance = read.table(file="mat_abundance_ab-sorensen.csv",sep=",")
simka_jaccard_abundance = read.table(file="mat_abundance_simka-jaccard.csv",sep=",")

chord_prevalence = read.table(file="mat_presenceAbsence_chord.csv",sep=",")
jaccard_prevalence = read.table(file="mat_presenceAbsence_jaccard.csv",sep=",")
kulczynski_prevalence = read.table(file="mat_presenceAbsence_kulczynski.csv",sep=",")
ochiai_prevalence = read.table(file="mat_presenceAbsence_ochiai.csv",sep=",")
whittaker_prevalence = read.table(file="mat_presenceAbsence_whittaker.csv",sep=",")
simka_jaccard_prevalence = read.table(file="mat_presenceAbsence_simka-jaccard.csv",sep=",")

jaccard_abundance = as.matrix(jaccard_abundance[-67,-67])
ab_jaccard_abundance = as.matrix(ab_jaccard_abundance[-67,-67])
braycurtis_abundance = as.matrix(braycurtis_abundance[-67,-67])
ab_ochiai_abundance = as.matrix(ab_ochiai_abundance[-67,-67])
ab_sorensen_abundance = as.matrix(ab_sorensen_abundance[-67,-67])
simka_jaccard_abundance = as.matrix(simka_jaccard_abundance[-67,-67])


chord_prevalence = as.matrix(chord_prevalence[-67,-67])
jaccard_prevalence = as.matrix(jaccard_prevalence[-67,-67])
kulczynski_prevalence = as.matrix(kulczynski_prevalence[-67,-67])
ochiai_prevalence = as.matrix(ochiai_prevalence[-67,-67])
whittaker_prevalence = as.matrix(whittaker_prevalence[-67,-67])
simka_jaccard_prevalence = as.matrix(simka_jaccard_prevalence[-67,-67])

k = 35
l = 5
brouillon = matrix(0,6,6)
# (k,l) = (30,2), (6,3), (40,4), (40,5)

# kNN1 = make.kNNG(jaccard_abundance, k = k, symm = TRUE, weight = FALSE)
# 
# kNN2 = make.kNNG(ab_jaccard_abundance, k = k, symm = TRUE, weight = FALSE)
# 
# kNN3 = make.kNNG(braycurtis_abundance, k = k, symm = TRUE, weight = FALSE)
# 
# kNN4 = make.kNNG(ab_ochiai_abundance, k = k, symm = TRUE, weight = FALSE)
# 
# kNN5 = make.kNNG(ab_sorensen_abundance, k = k, symm = TRUE, weight = FALSE)
# 
# kNN6 = make.kNNG(simka_jaccard_abundance, k = k, symm = TRUE, weight = FALSE)

kNN7 = make.kNNG(chord_prevalence, k = k, symm = TRUE, weight = FALSE)

kNN8 = make.kNNG(jaccard_prevalence, k = k, symm = TRUE, weight = FALSE)

kNN9 = make.kNNG(kulczynski_prevalence, k = k, symm = TRUE, weight = FALSE)

kNN10 = make.kNNG(ochiai_prevalence, k = k, symm = TRUE, weight = FALSE)

kNN11 = make.kNNG(whittaker_prevalence, k = k, symm = TRUE, weight = FALSE)

kNN12 = make.kNNG(simka_jaccard_prevalence, k = k, symm = TRUE, weight = FALSE)

index_one_vs_one = rep(0,100)
index_one_vs_all = rep(0,100)

for(nb_run in 1:100)
{
# res1 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# 
# res2 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# 
# res3 = spectral.clustering(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# 
# res4 = spectral.clustering(kNN4, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# 
# res5 = spectral.clustering(kNN5, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# 
# res6 = spectral.clustering(kNN6, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

res7 = spectral.clustering(kNN7, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

res8 = spectral.clustering(kNN8, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

res9 = spectral.clustering(kNN9, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

res10 = spectral.clustering(kNN10, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

res11 = spectral.clustering(kNN11, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

res12 = spectral.clustering(kNN12, normalised = TRUE, score = FALSE, K = l, adj = FALSE)


cluster = matrix(NA,nrow=210,ncol=12*l)
for(i in 1:l){
  # cluster[1:length(which(res1==i)),i] = which(res1==i)
  # cluster[1:length(which(res2==i)),i+l] = which(res2==i)
  # cluster[1:length(which(res3==i)),i+2*l] = which(res3==i)
  # cluster[1:length(which(res4==i)),i+3*l] = which(res4==i)
  # cluster[1:length(which(res5==i)),i+4*l] = which(res5==i)
  # cluster[1:length(which(res6==i)),i+5*l] = which(res6==i)
  cluster[1:length(which(res7==i)),i+6*l] = which(res7==i)
  cluster[1:length(which(res8==i)),i+7*l] = which(res8==i)
  cluster[1:length(which(res9==i)),i+8*l] = which(res9==i)
  cluster[1:length(which(res10==i)),i+9*l] = which(res10==i)
  cluster[1:length(which(res11==i)),i+10*l] = which(res11==i)
  cluster[1:length(which(res12==i)),i+11*l] = which(res12==i)
}

for(j in 7:11){
  index_one_vs_all[nb_run] = index_one_vs_all[nb_run] + adjustedRand(get(paste("res",7,sep="")),get(paste("res",j+1,sep="")),randMethod="Rand")
  for(m in (j+1):12){ 
    index_one_vs_one[nb_run] = index_one_vs_one[nb_run] + adjustedRand(get(paste("res",j,sep="")),get(paste("res",m,sep="")),randMethod="Rand")
  }
}

for(j in 7:12){
  for(m in 7:12){ 
    brouillon[j-6,m-6] = brouillon[j-6,m-6] + rand.index(get(paste("res",j,sep="")),get(paste("res",m,sep="")))/100
  }
}

index_one_vs_one[nb_run] = index_one_vs_one[nb_run]/15
index_one_vs_all[nb_run] = index_one_vs_all[nb_run]/5
cat(sprintf("Itération %s \n",nb_run))
}

mean(index_one_vs_one)
sd(index_one_vs_one)
mean(index_one_vs_all)
sd(index_one_vs_all)


