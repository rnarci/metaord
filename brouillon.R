rm(list=objects())

library(loe)
library(fcd)
library(clues)
library(cccd)
library(fossil)

# #Create a toy data
# x <- seq(-5,5,by=1)
# y <- seq(1,6,by=1)
# hx1 <- seq(-3.5,-1.5,by=0.5)
# hx2 <- seq(1.5,3.5,by=0.5)
# hy <- seq(2.5,4.5,by=0.5)
# D1 <- matrix(0,66,2)
# for(i in 1:11){
#   for(j in 1:6){
#     D1[i+11*(j-1),] <- c(x[i],y[j])
#   }
# }
# D2n <- matrix(0,25,2)
# D2p <- matrix(0,25,2)
# for(i in 1:5){
#   for(j in 1:5){
#     D2n[i+5*(j-1),] <- c(hx1[i],hy[j])
#     D2p[i+5*(j-1),] <- c(hx2[i],hy[j])
#   }
# }
# D2n <- D2n[-c(7,9,17,19),]
# D2p <- D2p[-c(7,9,17,19),]
# Data <- rbind(D1,D2n,D2p)
# Data <- scale(Data[order(Data[,1]),], scale=FALSE)
# 
# #Visualization of Data
# plot(Data,pch=20,xlab="",ylab="",cex=1,col=rainbow(108,start=.7,end=.9),
#      xlim=c(-7,7),ylim=c(-7,7))
# 
# #Creating a k-NN graph based on Data
# DM <- as.matrix(dist(Data))
# ADM <- make.kNNG(DM,k=5,symm=TRUE)
# 
# n_data = dim(DM)[1]
# 
# ################################################
# 
# plot(Data,pch=20,xlab="",ylab="",cex=1,col=rainbow(108,start=.7,end=.9),
#      xlim=c(-7,7),ylim=c(-7,7))
# 
# for(i in 1:(n_data-1)){
#   for(j in (i+1):n_data){
#     if(ADM[i,j]==1){
#       lines(c(Data[i,1],Data[j,1]),c(Data[i,2],Data[j,2]),col="black")
#       }
#   }
# }
# 
# res = spectral.clustering(ADM, normalised = TRUE, score = FALSE, K = 4, adj = FALSE)
# 
# plot(Data,xlab="",ylab="",cex=1,
#      xlim=c(-7,7),ylim=c(-7,7))
# 
# points(Data[which(res==1),1],Data[which(res==1),2],col="green")
# points(Data[which(res==2),1],Data[which(res==2),2],col="red")
# points(Data[which(res==3),1],Data[which(res==3),2],col="blue")
# points(Data[which(res==4),1],Data[which(res==4),2],col="black")
# 
# ####################################################
# 
# S = collect_all_statements_distance(DM)
# 
# ########## Optional : collect an arbitrary collection of statements from all statements
# 
# n_statements = 10000
# u = ceiling(runif(n_statements,0,dim(S)[1]))
# S_subset = S[u,]
# 
# ############################## Perform Algorithm 5 Clustering
# 
# k = 5
# l = 4
# sigma = 5
# 
# kRNG = k_RNG_clustering_unweighted(S=S_subset,n_data=n_data,k=k,l=l) 
# 
# # for(i in 1:(n_data-1)){
# #   for(j in (i+1):n_data){
# #     if(kRNG[i,j]==1){
# #       lines(c(data[i,1],data[j,1]),c(data[i,2],data[j,2]),col="blue")
# #       }
# #   }
# # }
# 
# res = spectral.clustering(kRNG, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# 
# plot(Data,xlab="",ylab="",cex=1,
#      xlim=c(-7,7),ylim=c(-7,7))
# 
# points(Data[which(res==1),1],Data[which(res==1),2],col="green")
# points(Data[which(res==2),1],Data[which(res==2),2],col="red")
# points(Data[which(res==3),1],Data[which(res==3),2],col="blue")
# points(Data[which(res==4),1],Data[which(res==4),2],col="black")







k = 2
l = 3
brouillon = matrix(0,6,6)
# fit1 = cmdscale(jaccard_abundance, eig=TRUE, k=k)
# 
# fit2 = cmdscale(ab_jaccard_abundance,eig=TRUE, k=k)
# 
# fit3 = cmdscale(braycurtis_abundance,eig=TRUE, k=k)
# 
# fit4 = cmdscale(ab_ochiai_abundance,eig=TRUE, k=k)
# 
# fit5 = cmdscale(ab_sorensen_abundance,eig=TRUE, k=k)
# 
# fit6 = cmdscale(simka_jaccard_abundance,eig=TRUE, k=k)

fit7 = cmdscale(chord_prevalence,eig=TRUE, k=k)

fit8 = cmdscale(jaccard_prevalence,eig=TRUE, k=k)

fit9 = cmdscale(kulczynski_prevalence,eig=TRUE, k=k)

fit10 = cmdscale(ochiai_prevalence,eig=TRUE, k=k)

fit11 = cmdscale(whittaker_prevalence,eig=TRUE, k=k)

fit12 = cmdscale(simka_jaccard_prevalence,eig=TRUE, k=k)


index_one_vs_one = rep(0,100)
index_one_vs_all = rep(0,100)

for(nb_run in 1:100)
{
# res1 = kmeans(fit1$points,l)$cluster
# 
# res2 = kmeans(fit2$points,l)$cluster
# 
# res3 = kmeans(fit3$points,l)$cluster
# 
# res4 = kmeans(fit4$points,l)$cluster
# 
# res5 = kmeans(fit5$points,l)$cluster
# 
# res6 = kmeans(fit6$points,l)$cluster

res7 = kmeans(fit7$points,l)$cluster

res8 = kmeans(fit8$points,l)$cluster

res9 = kmeans(fit9$points,l)$cluster

res10 = kmeans(fit10$points,l)$cluster

res11 = kmeans(fit11$points,l)$cluster

res12 = kmeans(fit12$points,l)$cluster


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













































k = 30
l = 2
sigma = 0.5
brouillon = matrix(0,12,12)
kRNG1 = make.kRNG(jaccard_abundance, k = k)

kRNG2 = make.kRNG(ab_jaccard_abundance, k = k)

kRNG3 = make.kRNG(braycurtis_abundance, k = k)

kRNG4 = make.kRNG(ab_ochiai_abundance, k = k)

kRNG5 = make.kRNG(ab_sorensen_abundance, k = k)

kRNG6 = make.kRNG(simka_jaccard_abundance, k = k)

kRNG7 = make.kRNG(chord_prevalence, k = k)

kRNG8 = make.kRNG(jaccard_prevalence, k = k)

kRNG9 = make.kRNG(kulczynski_prevalence, k = k)

kRNG10 = make.kRNG(ochiai_prevalence, k = k)

kRNG11 = make.kRNG(whittaker_prevalence, k = k)

kRNG12 = make.kRNG(simka_jaccard_prevalence, k = k)

index_one_vs_one = rep(0,20)
index_one_vs_all = rep(0,20)

for(nb_run in 1:20)
{
  res1 = spectral.clustering(kRNG1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

  res2 = spectral.clustering(kRNG2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

  res3 = spectral.clustering(kRNG3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

  res4 = spectral.clustering(kRNG4, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

  res5 = spectral.clustering(kRNG5, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

  res6 = spectral.clustering(kRNG6, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  res7 = spectral.clustering(kRNG7, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  res8 = spectral.clustering(kRNG8, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  res9 = spectral.clustering(kRNG9, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  res10 = spectral.clustering(kRNG10, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  res11 = spectral.clustering(kRNG11, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  res12 = spectral.clustering(kRNG12, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  
  cluster = matrix(NA,nrow=210,ncol=12*l)
  for(i in 1:l){
    cluster[1:length(which(res1==i)),i] = which(res1==i)
    cluster[1:length(which(res2==i)),i+l] = which(res2==i)
    cluster[1:length(which(res3==i)),i+2*l] = which(res3==i)
    cluster[1:length(which(res4==i)),i+3*l] = which(res4==i)
    cluster[1:length(which(res5==i)),i+4*l] = which(res5==i)
    cluster[1:length(which(res6==i)),i+5*l] = which(res6==i)
    cluster[1:length(which(res7==i)),i+6*l] = which(res7==i)
    cluster[1:length(which(res8==i)),i+7*l] = which(res8==i)
    cluster[1:length(which(res9==i)),i+8*l] = which(res9==i)
    cluster[1:length(which(res10==i)),i+9*l] = which(res10==i)
    cluster[1:length(which(res11==i)),i+10*l] = which(res11==i)
    cluster[1:length(which(res12==i)),i+11*l] = which(res12==i)
  }
  
  for(j in 1:11){
    index_one_vs_all[nb_run] = index_one_vs_all[nb_run] + adjustedRand(get(paste("res",1,sep="")),get(paste("res",j+1,sep="")),randMethod="Rand")
    for(m in (j+1):12){ 
      index_one_vs_one[nb_run] = index_one_vs_one[nb_run] + adjustedRand(get(paste("res",j,sep="")),get(paste("res",m,sep="")),randMethod="Rand")
    }
  }
  for(j in 1:12){
    for(m in 1:12){ 
      brouillon[j,m] = brouillon[j,m] + rand.index(get(paste("res",j,sep="")),get(paste("res",m,sep="")))/20
    }
  }
  index_one_vs_one[nb_run] = index_one_vs_one[nb_run]/66 
  index_one_vs_all[nb_run] = index_one_vs_all[nb_run]/11 
  cat(sprintf("Itération %s \n",nb_run))
}

mean(index_one_vs_one)
sd(index_one_vs_one)
mean(index_one_vs_all)
sd(index_one_vs_all)


















k = 7
l = 9

fit1 = isoMDS(jaccard_abundance, k=k)

fit2 = isoMDS(ab_jaccard_abundance, k=k)

fit3 = isoMDS(braycurtis_abundance, k=k)

fit4 = isoMDS(ab_ochiai_abundance, k=k)

fit5 = isoMDS(ab_sorensen_abundance, k=k)

fit6 = isoMDS(simka_jaccard_abundance, k=k)

fit7 = isoMDS(chord_prevalence, k=k)

fit8 = isoMDS(jaccard_prevalence, k=k)

fit9 = isoMDS(kulczynski_prevalence, k=k)

fit10 = isoMDS(ochiai_prevalence, k=k)

fit11 = isoMDS(whittaker_prevalence, k=k)

fit12 = isoMDS(simka_jaccard_prevalence, k=k)


index_one_vs_one = rep(0,20)
index_one_vs_all = rep(0,20)

for(nb_run in 1:20)
{
  res1 = kmeans(fit1$points,l)$cluster

  res2 = kmeans(fit2$points,l)$cluster

  res3 = kmeans(fit3$points,l)$cluster

  res4 = kmeans(fit4$points,l)$cluster

  res5 = kmeans(fit5$points,l)$cluster

  res6 = kmeans(fit6$points,l)$cluster
  
  res7 = kmeans(fit7$points,l)$cluster
  
  res8 = kmeans(fit8$points,l)$cluster
  
  res9 = kmeans(fit9$points,l)$cluster
  
  res10 = kmeans(fit10$points,l)$cluster
  
  res11 = kmeans(fit11$points,l)$cluster
  
  res12 = kmeans(fit12$points,l)$cluster
  
  
  cluster = matrix(NA,nrow=210,ncol=12*l)
  for(i in 1:l){
    cluster[1:length(which(res1==i)),i] = which(res1==i)
    cluster[1:length(which(res2==i)),i+l] = which(res2==i)
    cluster[1:length(which(res3==i)),i+2*l] = which(res3==i)
    cluster[1:length(which(res4==i)),i+3*l] = which(res4==i)
    cluster[1:length(which(res5==i)),i+4*l] = which(res5==i)
    cluster[1:length(which(res6==i)),i+5*l] = which(res6==i)
    cluster[1:length(which(res7==i)),i+6*l] = which(res7==i)
    cluster[1:length(which(res8==i)),i+7*l] = which(res8==i)
    cluster[1:length(which(res9==i)),i+8*l] = which(res9==i)
    cluster[1:length(which(res10==i)),i+9*l] = which(res10==i)
    cluster[1:length(which(res11==i)),i+10*l] = which(res11==i)
    cluster[1:length(which(res12==i)),i+11*l] = which(res12==i)
  }
  
  for(j in 1:11){
    index_one_vs_all[nb_run] = index_one_vs_all[nb_run] + adjustedRand(get(paste("res",1,sep="")),get(paste("res",j+1,sep="")),randMethod="Rand")
    for(m in (j+1):12){
      index_one_vs_one[nb_run] = index_one_vs_one[nb_run] + adjustedRand(get(paste("res",j,sep="")),get(paste("res",m,sep="")),randMethod="Rand")
    }
  }
  for(j in 1:12){
    for(m in 1:12){ 
      brouillon[j,m] = brouillon[j,m] + adjustedRand(get(paste("res",j,sep="")),get(paste("res",m,sep="")),randMethod="Rand")/20
    }
  }
  index_one_vs_one[nb_run] = index_one_vs_one[nb_run]/66
  index_one_vs_all[nb_run] = index_one_vs_all[nb_run]/11
  cat(sprintf("Itération %s \n",nb_run))
}

mean(index_one_vs_one)
sd(index_one_vs_one)
mean(index_one_vs_all)
sd(index_one_vs_all)











































k = 20
l = 9
brouillon = matrix(0,12,12)
# (k,l) = (30,2), (6,3), (40,4), (40,5)

kNN1 = make.kNNG(jaccard_abundance, k = k, symm = TRUE, weight = TRUE)

kNN2 = make.kNNG(ab_jaccard_abundance, k = k, symm = TRUE, weight = TRUE)

kNN3 = make.kNNG(braycurtis_abundance, k = k, symm = TRUE, weight = TRUE)

kNN4 = make.kNNG(ab_ochiai_abundance, k = k, symm = TRUE, weight = TRUE)

kNN5 = make.kNNG(ab_sorensen_abundance, k = k, symm = TRUE, weight = TRUE)

kNN6 = make.kNNG(simka_jaccard_abundance, k = k, symm = TRUE, weight = TRUE)

kNN7 = make.kNNG(chord_prevalence, k = k, symm = TRUE, weight = TRUE)

kNN8 = make.kNNG(jaccard_prevalence, k = k, symm = TRUE, weight = TRUE)

kNN9 = make.kNNG(kulczynski_prevalence, k = k, symm = TRUE, weight = TRUE)

kNN10 = make.kNNG(ochiai_prevalence, k = k, symm = TRUE, weight = TRUE)

kNN11 = make.kNNG(whittaker_prevalence, k = k, symm = TRUE, weight = TRUE)

kNN12 = make.kNNG(simka_jaccard_prevalence, k = k, symm = TRUE, weight = TRUE)

index_one_vs_one = rep(0,20)
index_one_vs_all = rep(0,20)

for(nb_run in 1:20)
{
  res1 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

  res2 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

  res3 = spectral.clustering(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

  res4 = spectral.clustering(kNN4, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

  res5 = spectral.clustering(kNN5, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

  res6 = spectral.clustering(kNN6, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  res7 = spectral.clustering(kNN7, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  res8 = spectral.clustering(kNN8, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  res9 = spectral.clustering(kNN9, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  res10 = spectral.clustering(kNN10, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  res11 = spectral.clustering(kNN11, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  res12 = spectral.clustering(kNN12, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  
  
  cluster = matrix(NA,nrow=210,ncol=12*l)
  for(i in 1:l){
    cluster[1:length(which(res1==i)),i] = which(res1==i)
    cluster[1:length(which(res2==i)),i+l] = which(res2==i)
    cluster[1:length(which(res3==i)),i+2*l] = which(res3==i)
    cluster[1:length(which(res4==i)),i+3*l] = which(res4==i)
    cluster[1:length(which(res5==i)),i+4*l] = which(res5==i)
    cluster[1:length(which(res6==i)),i+5*l] = which(res6==i)
    cluster[1:length(which(res7==i)),i+6*l] = which(res7==i)
    cluster[1:length(which(res8==i)),i+7*l] = which(res8==i)
    cluster[1:length(which(res9==i)),i+8*l] = which(res9==i)
    cluster[1:length(which(res10==i)),i+9*l] = which(res10==i)
    cluster[1:length(which(res11==i)),i+10*l] = which(res11==i)
    cluster[1:length(which(res12==i)),i+11*l] = which(res12==i)
  }
  
  for(j in 1:11){
    index_one_vs_all[nb_run] = index_one_vs_all[nb_run] + adjustedRand(get(paste("res",1,sep="")),get(paste("res",j+1,sep="")),randMethod="Rand")
    for(m in (j+1):12){ 
      index_one_vs_one[nb_run] = index_one_vs_one[nb_run] + adjustedRand(get(paste("res",j,sep="")),get(paste("res",m,sep="")),randMethod="Rand")
    }
  }
  for(j in 1:12){
    for(m in 1:12){ 
     brouillon[j,m] = brouillon[j,m] + adjustedRand(get(paste("res",j,sep="")),get(paste("res",m,sep="")),randMethod="Rand")/20
    }
  }
  index_one_vs_one[nb_run] = index_one_vs_one[nb_run]/66 
  index_one_vs_all[nb_run] = index_one_vs_all[nb_run]/11
  cat(sprintf("Itération %s \n",nb_run))
}

mean(index_one_vs_one)
sd(index_one_vs_one)
mean(index_one_vs_all)
sd(index_one_vs_all)





























l=5










kNN11 = make.kNNG(chord_prevalence, k = k, symm = TRUE, weight = FALSE)
kNN8 = make.kNNG(jaccard_prevalence, k = k, symm = TRUE, weight = FALSE)
res11 = spectral.clustering(kNN11, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res8 = spectral.clustering(kNN8, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res11 = specc(kNN11, centers=l)[1:140]
# res8 = specc(kNN8,centers=l)[1:140]
adjustedRand(get(paste("res",11,sep="")),get(paste("res",8,sep="")),randMethod="Rand")

# kNN11 = make.kNNG(chord_prevalence, k = 35, symm = TRUE, weight = FALSE)
# kNN8 = make.kNNG(jaccard_prevalence, k = 35, symm = TRUE, weight = FALSE)
# kRNG11 = make.kRNG(chord_prevalence, k = k)
# kRNG8 = make.kRNG(jaccard_prevalence, k = k)
# res11 = spectral.clustering(kRNG11, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res8 = spectral.clustering(kRNG8, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# adjustedRand(get(paste("res",11,sep="")),get(paste("res",8,sep="")),randMethod="Rand")

resbis11 = spectral.clustering.rw(kNN11, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
resbis8 = spectral.clustering.rw(kNN8, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
adjustedRand(get(paste("resbis",11,sep="")),get(paste("resbis",8,sep="")),randMethod="Rand")

fit11 = cmdscale(chord_prevalence,eig=TRUE, k=4)
fit8 = cmdscale(jaccard_prevalence,eig=TRUE, k=4)
resbisbis11 = kmeans(fit11$points,l)$cluster
resbisbis8 = kmeans(fit8$points,l)$cluster
adjustedRand(get(paste("resbisbis",11,sep="")),get(paste("resbisbis",8,sep="")),randMethod="Rand")


adjustedRand(get(paste("resbis",11,sep="")),get(paste("res",11,sep="")),randMethod="Rand")
adjustedRand(get(paste("res",11,sep="")),get(paste("resbisbis",11,sep="")),randMethod="Rand")
adjustedRand(get(paste("resbis",11,sep="")),get(paste("resbisbis",11,sep="")),randMethod="Rand")

cat(sprintf("\n ARI spectral.clustering = %s \n ARI spectral.clustering.rw = %s \n ARI kmeans = %s \n ARI between spectral.clustering and spectral.clustering.rw = %s \n ARI between spectral.clustering and kmeans = %s \n ARI between spectral.clustering.rw and kmeans = %s",
    round(adjustedRand(get(paste("res",11,sep="")),get(paste("res",8,sep="")),randMethod="Rand"),3), round(adjustedRand(get(paste("resbis",11,sep="")),get(paste("resbis",8,sep="")),randMethod="Rand"),3),
          round(adjustedRand(get(paste("resbisbis",11,sep="")),get(paste("resbisbis",8,sep="")),randMethod="Rand"),3),round(adjustedRand(get(paste("resbis",11,sep="")),get(paste("res",11,sep="")),randMethod="Rand"),3),
                round(adjustedRand(get(paste("res",11,sep="")),get(paste("resbisbis",11,sep="")),randMethod="Rand"),3),round(adjustedRand(get(paste("resbis",11,sep="")),get(paste("resbisbis",11,sep="")),randMethod="Rand"),3)))






kNN11 = make.kNNG(chord_prevalence, k = k, symm = TRUE, weight = FALSE)
ktildeNN11 = make.kNNG(chord_prevalence, k = ktilde, symm = TRUE, weight = FALSE)
res11 = specc(kNN11, centers=l)[1:140]
restilde11 = specc(ktildeNN11,centers=l)[1:140]
adjustedRand(get(paste("res",11,sep="")),get(paste("restilde",11,sep="")),randMethod="Rand")

kNN11 = make.kNNG(chord_prevalence, k = k, symm = TRUE, weight = FALSE)
ktildeNN11 = make.kNNG(chord_prevalence, k = ktilde, symm = TRUE, weight = FALSE)
res11 = spectral.clustering.rw(kNN11, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
restilde11 = spectral.clustering.rw(ktildeNN11, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
adjustedRand(get(paste("res",11,sep="")),get(paste("restilde",11,sep="")),randMethod="Rand")















































count_neighbors=c()
num_matrix = 4
for(i in 1:n_data){
    count_neighbors[i] = sum(kNN[[num_matrix]][i,]) + sum(kNN[[num_matrix]][,i])
    }
count_neighbors-k

