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


















































rm(list=objects())

library(psych)

Y = matrix(c(1,2,4,7,8,9,3,4,2,5,6,4,1,5,8,9,5,7,1,2),4,5)
Y = scale(Y)
Z=Y%*%t(Y)
D = as.matrix(dist(Y,method="euclidean",diag=TRUE))

A = -0.5*D^2

Id = diag(identity(4))
un = t(t(c(1,1,1,1)))
scale = Id-(1/4)*un%*%t(un)

G = scale%*%A%*%scale


































########## Au cas ou


# rm(list=objects())
# 
# library(fcd)
# library(loe)
# library(clues)
# library(cccd)
# library(fossil)
# 
# jaccard_abundance = read.table(file="mat_abundance_jaccard.csv",sep=",")
# ab_jaccard_abundance = read.table(file="mat_abundance_ab-jaccard.csv",sep=",")
# braycurtis_abundance = read.table(file="mat_abundance_braycurtis.csv",sep=",")
# ab_ochiai_abundance = read.table(file="mat_abundance_ab-ochiai.csv",sep=",")
# ab_sorensen_abundance = read.table(file="mat_abundance_ab-sorensen.csv",sep=",")
# simka_jaccard_abundance = read.table(file="mat_abundance_simka-jaccard.csv",sep=",")
# 
# chord_prevalence = read.table(file="mat_presenceAbsence_chord.csv",sep=",")
# jaccard_prevalence = read.table(file="mat_presenceAbsence_jaccard.csv",sep=",")
# kulczynski_prevalence = read.table(file="mat_presenceAbsence_kulczynski.csv",sep=",")
# ochiai_prevalence = read.table(file="mat_presenceAbsence_ochiai.csv",sep=",")
# whittaker_prevalence = read.table(file="mat_presenceAbsence_whittaker.csv",sep=",")
# simka_jaccard_prevalence = read.table(file="mat_presenceAbsence_simka-jaccard.csv",sep=",")
# 
# jaccard_abundance = as.matrix(jaccard_abundance[-67,-67])
# ab_jaccard_abundance = as.matrix(ab_jaccard_abundance[-67,-67])
# braycurtis_abundance = as.matrix(braycurtis_abundance[-67,-67])
# ab_ochiai_abundance = as.matrix(ab_ochiai_abundance[-67,-67])
# ab_sorensen_abundance = as.matrix(ab_sorensen_abundance[-67,-67])
# simka_jaccard_abundance = as.matrix(simka_jaccard_abundance[-67,-67])
# 
# 
# chord_prevalence = as.matrix(chord_prevalence[-67,-67])
# jaccard_prevalence = as.matrix(jaccard_prevalence[-67,-67])
# kulczynski_prevalence = as.matrix(kulczynski_prevalence[-67,-67])
# ochiai_prevalence = as.matrix(ochiai_prevalence[-67,-67])
# whittaker_prevalence = as.matrix(whittaker_prevalence[-67,-67])
# simka_jaccard_prevalence = as.matrix(simka_jaccard_prevalence[-67,-67])
# 
# k = 35
# l = 2
# brouillon = matrix(0,6,6)
# # (k,l) = (30,2), (6,3), (40,4), (40,5)
# 
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
# 
# kNN7 = make.kNNG(chord_prevalence, k = k, symm = TRUE, weight = FALSE)
# 
# kNN8 = make.kNNG(jaccard_prevalence, k = k, symm = TRUE, weight = FALSE)
# 
# kNN9 = make.kNNG(kulczynski_prevalence, k = k, symm = TRUE, weight = FALSE)
# 
# kNN10 = make.kNNG(ochiai_prevalence, k = k, symm = TRUE, weight = FALSE)
# 
# kNN11 = make.kNNG(whittaker_prevalence, k = k, symm = TRUE, weight = FALSE)
# 
# kNN12 = make.kNNG(simka_jaccard_prevalence, k = k, symm = TRUE, weight = FALSE)
# 
# index_one_vs_one = rep(0,100)
# index_one_vs_all = rep(0,100)
# 
# for(nb_run in 1:100)
# {
#   res1 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   
#   res2 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   
#   res3 = spectral.clustering(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   
#   res4 = spectral.clustering(kNN4, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   
#   res5 = spectral.clustering(kNN5, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   
#   res6 = spectral.clustering(kNN6, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   
#   res7 = spectral.clustering(kNN7, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   
#   res8 = spectral.clustering(kNN8, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   
#   res9 = spectral.clustering(kNN9, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   
#   res10 = spectral.clustering(kNN10, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   
#   res11 = spectral.clustering(kNN11, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   
#   res12 = spectral.clustering(kNN12, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   
#   
#   cluster = matrix(NA,nrow=210,ncol=12*l)
#   for(i in 1:l){
#     cluster[1:length(which(res1==i)),i] = which(res1==i)
#     cluster[1:length(which(res2==i)),i+l] = which(res2==i)
#     cluster[1:length(which(res3==i)),i+2*l] = which(res3==i)
#     cluster[1:length(which(res4==i)),i+3*l] = which(res4==i)
#     cluster[1:length(which(res5==i)),i+4*l] = which(res5==i)
#     cluster[1:length(which(res6==i)),i+5*l] = which(res6==i)
#     cluster[1:length(which(res7==i)),i+6*l] = which(res7==i)
#     cluster[1:length(which(res8==i)),i+7*l] = which(res8==i)
#     cluster[1:length(which(res9==i)),i+8*l] = which(res9==i)
#     cluster[1:length(which(res10==i)),i+9*l] = which(res10==i)
#     cluster[1:length(which(res11==i)),i+10*l] = which(res11==i)
#     cluster[1:length(which(res12==i)),i+11*l] = which(res12==i)
#   }
#   
#   for(j in 7:11){
#     index_one_vs_all[nb_run] = index_one_vs_all[nb_run] + adjustedRand(get(paste("res",7,sep="")),get(paste("res",j+1,sep="")),randMethod="Rand")
#     for(m in (j+1):12){ 
#       index_one_vs_one[nb_run] = index_one_vs_one[nb_run] + adjustedRand(get(paste("res",j,sep="")),get(paste("res",m,sep="")),randMethod="Rand")
#     }
#   }
#   
#   for(j in 7:12){
#     for(m in 7:12){ 
#       brouillon[j-6,m-6] = brouillon[j-6,m-6] + rand.index(get(paste("res",j,sep="")),get(paste("res",m,sep="")))/100
#     }
#   }
#   
#   index_one_vs_one[nb_run] = index_one_vs_one[nb_run]/15
#   index_one_vs_all[nb_run] = index_one_vs_all[nb_run]/5
#   cat(sprintf("Itération %s \n",nb_run))
# }
# 
# mean(index_one_vs_one)
# sd(index_one_vs_one)
# mean(index_one_vs_all)
# sd(index_one_vs_all)
































############################## Choose size fraction

size_fraction = "0.8-5"

if(size_fraction=="0.8-5"){
  setwd(dir="/home/rnarci/Bureau/CDDrnarci/Donnees/simka_matrices_17-05-22/Simka_0.8-5") 
}

if(size_fraction=="5-20"){
  setwd(dir="/home/rnarci/Bureau/CDDrnarci/Donnees/simka_matrices_17-05-22/Simka_5-20") 
}

if(size_fraction=="180-2000"){
  setwd(dir="/home/rnarci/Bureau/CDDrnarci/Donnees/simka_matrices_17-05-22/Simka_180-2000") 
}

if(size_fraction=="0.22-3"){
  setwd(dir="/home/rnarci/Bureau/CDDrnarci/Donnees/simka_matrices_17-05-22/Simka_0.22-3") 
}

if(size_fraction=="0-0.2"){
  setwd(dir="/home/rnarci/Bureau/CDDrnarci/Donnees/simka_matrices_17-05-22/Simka_0-0.2") 
}

if(size_fraction=="20-180"){
  setwd(dir="/home/rnarci/Bureau/CDDrnarci/Donnees/simka_matrices_17-05-22/Simka_20-180") 
}

############################## Import data (note : pour l'instant, je ne prends pas les matrices de distance asymetriques et celle de sorensen-braycurtis PA)

if(size_fraction=="0.8-5" | size_fraction=="5-20" | size_fraction=="180-2000"){
  jaccard_abundance2 = read.table(file="mat_abundance_jaccard.csv",sep="")
  if(dim(jaccard_abundance2)[2]==1){
    jaccard_abundance2 = read.table(file="mat_abundance_jaccard.csv",sep=";")
    sep = ";"
  }else{
    sep = ""
  }
  delete=c(which(duplicated(jaccard_abundance2)==TRUE)) # remove duplicates
  jaccard_abundance = jaccard_abundance2[-delete,-delete]
  sample = jaccard_abundance2[-1,1]
  
  jaccard_abundance2 = read.table(file="mat_abundance_jaccard.csv",sep=sep)
  ab_jaccard_abundance2 = read.table(file="mat_abundance_ab-jaccard.csv",sep=sep)
  braycurtis_abundance2 = read.table(file="mat_abundance_braycurtis.csv",sep=sep)
  ab_ochiai_abundance2 = read.table(file="mat_abundance_ab-ochiai.csv",sep=sep)
  ab_sorensen_abundance2 = read.table(file="mat_abundance_ab-sorensen.csv",sep=sep)
  simka_jaccard_abundance2 = read.table(file="mat_abundance_simka-jaccard.csv",sep=sep)
  
  chord_prevalence2 = read.table(file="mat_presenceAbsence_chord.csv",sep=sep)
  jaccard_prevalence2 = read.table(file="mat_presenceAbsence_jaccard.csv",sep=sep)
  kulczynski_prevalence2 = read.table(file="mat_presenceAbsence_kulczynski.csv",sep=sep)
  ochiai_prevalence2 = read.table(file="mat_presenceAbsence_ochiai.csv",sep=sep)
  whittaker_prevalence2 = read.table(file="mat_presenceAbsence_whittaker.csv",sep=sep)
  simka_jaccard_prevalence2 = read.table(file="mat_presenceAbsence_simka-jaccard.csv",sep=sep)
  # sorensen_braycurtis_prevalence2 = read.table(file="mat_presenceAbsence_sorensen-braycurtis.csv",sep="") doesn't work
  
  ############### Distance matrices
  
  jaccard_abundance2 = unname(as.matrix(jaccard_abundance2)[c(-1,-delete),c(-1,-delete)])
  ab_jaccard_abundance2 = unname(as.matrix(ab_jaccard_abundance2)[c(-1,-delete),c(-1,-delete)])
  braycurtis_abundance2 = unname(as.matrix(braycurtis_abundance2)[c(-1,-delete),c(-1,-delete)])
  ab_ochiai_abundance2 = unname(as.matrix(ab_ochiai_abundance2)[c(-1,-delete),c(-1,-delete)])
  ab_sorensen_abundance2 = unname(as.matrix(ab_sorensen_abundance2)[c(-1,-delete),c(-1,-delete)])
  simka_jaccard_abundance2 = unname(as.matrix(simka_jaccard_abundance2)[c(-1,-delete),c(-1,-delete)])
  
  
  chord_prevalence2 = unname(as.matrix(chord_prevalence2)[c(-1,-delete),c(-1,-delete)])
  jaccard_prevalence2 = unname(as.matrix(jaccard_prevalence2)[c(-1,-delete),c(-1,-delete)])
  kulczynski_prevalence2 = unname(as.matrix(kulczynski_prevalence2)[c(-1,-delete),c(-1,-delete)])
  ochiai_prevalence2 = unname(as.matrix(ochiai_prevalence2)[c(-1,-delete),c(-1,-delete)])
  whittaker_prevalence2 = unname(as.matrix(whittaker_prevalence2)[c(-1,-delete),c(-1,-delete)])
  simka_jaccard_prevalence2 = unname(as.matrix(simka_jaccard_prevalence2)[c(-1,-delete),c(-1,-delete)])
  # sorensen_braycurtis_prevalence2 = unname(as.matrix(sorensen_braycurtis_prevalence2)[c(-1,-delete),c(-1,-delete)])
  
  n = dim(jaccard_abundance2)[1]
  
  jaccard_abundance = matrix(NA,nrow=n,ncol=n)
  ab_jaccard_abundance = matrix(NA,nrow=n,ncol=n)
  braycurtis_abundance = matrix(NA,nrow=n,ncol=n)
  ab_ochiai_abundance = matrix(NA,nrow=n,ncol=n)
  ab_sorensen_abundance = matrix(NA,nrow=n,ncol=n)
  simka_jaccard_abundance = matrix(NA,nrow=n,ncol=n)
  
  chord_prevalence = matrix(NA,nrow=n,ncol=n)
  jaccard_prevalence = matrix(NA,nrow=n,ncol=n)
  kulczynski_prevalence = matrix(NA,nrow=n,ncol=n)
  ochiai_prevalence = matrix(NA,nrow=n,ncol=n)
  whittaker_prevalence = matrix(NA,nrow=n,ncol=n)
  simka_jaccard_prevalence = matrix(NA,nrow=n,ncol=n)
  # sorensen_braycurtis_prevalence = matrix(NA,nrow=n,ncol=n)
  
  for(j in 1:n){
    jaccard_abundance[,j] = as.numeric(jaccard_abundance2[,j])
    ab_jaccard_abundance[,j] = as.numeric(ab_jaccard_abundance2[,j])
    braycurtis_abundance[,j] = as.numeric(braycurtis_abundance2[,j])
    ab_ochiai_abundance[,j] = as.numeric(ab_ochiai_abundance2[,j])
    ab_sorensen_abundance[,j] = as.numeric(ab_sorensen_abundance2[,j])
    simka_jaccard_abundance[,j] = as.numeric(simka_jaccard_abundance2[,j])
    
    chord_prevalence[,j] = as.numeric(chord_prevalence2[,j])
    jaccard_prevalence[,j] = as.numeric(jaccard_prevalence2[,j])
    kulczynski_prevalence[,j] = as.numeric(kulczynski_prevalence2[,j])
    ochiai_prevalence[,j] = as.numeric(ochiai_prevalence2[,j])
    whittaker_prevalence[,j] = as.numeric(whittaker_prevalence2[,j])
    simka_jaccard_prevalence[,j] = as.numeric(simka_jaccard_prevalence2[,j])
    # sorensen_braycurtis_prevalence[,j] = as.numeric(sorensen_braycurtis_prevalence2[,j])
  }
  
}else{
  sep = ""
  jaccard_abundance = read.table(file="mat_abundance_jaccard.csv",sep=sep)
  ochiai_abundance = read.table(file="mat_abundance_ochiai.csv",sep=sep)
  sorensen_abundance = read.table(file="mat_abundance_sorensen.csv",sep=sep)
  simka_jaccard_abundance = read.table(file="mat_abundance_simka-jaccard.csv",sep=sep)
  
  chord_hellinger_prevalence = read.table(file="mat_presenceAbsence_chord-hellinger.csv",sep=sep)
  jaccard_canberra_prevalence = read.table(file="mat_presenceAbsence_jaccard-canberra.csv",sep=sep)
  kulczynski_prevalence = read.table(file="mat_presenceAbsence_kulczynski.csv",sep=sep)
  ochiai_prevalence = read.table(file="mat_presenceAbsence_ochiai.csv",sep=sep)
  whittaker_prevalence = read.table(file="mat_presenceAbsence_whittaker.csv",sep=sep)
  simka_jaccard_prevalence = read.table(file="mat_presenceAbsence_simka-jaccard.csv",sep=sep)
  sorensen_braycurtis_prevalence = read.table(file="mat_presenceAbsence_sorensen-braycurtis.csv",sep=sep)
  
  ############### Distance matrices
  
  jaccard_abundance = unname(as.matrix(jaccard_abundance))
  ochiai_abundance = unname(as.matrix(ochiai_abundance))
  sorensen_abundance = unname(as.matrix(sorensen_abundance))
  simka_jaccard_abundance = unname(as.matrix(simka_jaccard_abundance))
  
  
  chord_hellinger_prevalence = unname(as.matrix(chord_hellinger_prevalence))
  jaccard_canberra_prevalence = unname(as.matrix(jaccard_canberra_prevalence))
  kulczynski_prevalence = unname(as.matrix(kulczynski_prevalence))
  ochiai_prevalence = unname(as.matrix(ochiai_prevalence))
  whittaker_prevalence = unname(as.matrix(whittaker_prevalence))
  simka_jaccard_prevalence = unname(as.matrix(simka_jaccard_prevalence))
  sorensen_braycurtis_prevalence = unname(as.matrix(sorensen_braycurtis_prevalence))
}



# ## Imputation de donnees manquantes
# 
# ## Methode 1 :
# summary(design)
# symnum(cor(design, use = "complete.obs"))
# 
# missDummy <- function(t)
# {
#   x <- dim(length(t)) 
#   x[which(!is.na(t))] = 1
#   x[which(is.na(t))] = 0
#   return(x)
# }
# 
# #Now we will use this function to create a dummy variable that will indicate missing value using 0, otherwise willtake the value 1.
# 
# design_bis = design
# design_bis$dummy <- missDummy(design$Phosphates)
# 
# #Let us take a look at the data
# design_bis
# 
# lm(Phosphates ~ Silicate, design)
# 
# for(i in 1:nrow(design))
# {
#   if(design_bis$dummy[i] == 0)
#   {
#     # design$Phosphates[i] = -1.184303 + 0.009963*design$Silicate[i] + 0.002018*design$SSD[i]
#     design$Phosphates[i] = 0.23558 + 0.02207*design$Silicate[i]
#   }
# }
# 
# design$Phosphates
# 
# ## Methode 2 :
# 
# # library(imputeR)
# library(mice)
# 
# # impute(design)
# complete(mice(design,method="norm.predict",m=1))$Phosphates

































######### Distances entre valeurs predites



rm(list=objects())

Y = matrix(c(1,2,4,7,8,9,3,4,2,5,6,4,1,5,8,9,5,7,1,2),4,5)
Y = scale(Y,scale=F)

X = matrix(c(4,5,3,1,7,8,4,9,1,2,5,6),4,3)
H = X%*%solve(t(X)%*%X)%*%t(X)

# Essai 1

v1 = sum(((H%*%Y)[1,])^2)

D = as.matrix(dist(Y,method="euclidean",diag=TRUE))
A = -0.5*D^2
n = 4
J = diag(rep(1,n)) - matrix(1,n,n)/n
G = J%*%A%*%J
v2 = H[1,]%*%G%*%H[,1]

# Essai 2 

Y_hat = H%*%Y 
R = Y - Y_hat
# D1 = matrix(NA,4,4)
# 
# for(i in 1:4){
#   for(j in 1:4){
#     D1[i,j] = sum((Y_hat[i,]-Y_hat[j,])^2)    
#   }
# }

D2 = matrix(NA,4,4)

for(i in 1:4){
  for(j in 1:4){
    D2[i,j] = H[i,]%*%G%*%H[,i] + H[j,]%*%G%*%H[,j] - 2*H[i,]%*%G%*%H[,j]
  }
}

I = diag(rep(1,n))
D1 = matrix(NA,4,4)

for(i in 1:4){
  for(j in 1:4){
    D1[i,j] = (I-H)[i,]%*%G%*%(I-H)[,i] + (I-H)[j,]%*%G%*%(I-H)[,j] - 2*(I-H)[i,]%*%G%*%(I-H)[,j]
  }
}

######## D_R en fonction de D_P, G et H

test = matrix(NA,4,4)
test1 = matrix(NA,4,4)
test2 = matrix(NA,4,4)

for(i in 1:4){
  for(j in 1:4){
    # test[i,j] = G[i,i] + G[j,j] -2*G[i,j] + D2[i,j] + 2*(H[i,]-H[j,])%*%(G[j,]-G[i,])
    test[i,j] = D[i,j]^2 + D2[i,j] + 2*(H[i,]-H[j,])%*%(G[j,]-G[i,])
  }
}

D_R = matrix(NA,4,4)
D_P = matrix(NA,4,4)

for(i in 1:4){
  for(j in 1:4){
    # test1[i,j] = D1[i,j] + D2[i,j]
    # test2[i,j] = D1[i,j] - D2[i,j] - 2*(H[i,]-H[j,])%*%(G[j,]-G[i,])
    D_R[i,j] = D1[i,j] 
    D_P[i,j] = -D2[i,j] - 2*(H[i,]-H[j,])%*%(G[j,]-G[i,])
  }
}
test = D_R + D_P # test = D^2


######### Distances entre résidus



rm(list=objects())

Y = matrix(c(1,2,4,7,8,9,3,4,2,5,6,4,1,5,8,9,5,7,1,2),4,5)
Y = scale(Y)

X = matrix(c(4,5,3,1,7,8,4,9,1,2,5,6),4,3)
H = X%*%solve(t(X)%*%X)%*%t(X)



D = as.matrix(dist(Y,method="euclidean",diag=TRUE))
A = -0.5*D^2
n = 4
J = diag(rep(1,n)) - matrix(1,n,n)/n
G = J%*%A%*%J

Y_hat = H%*%Y 
R = Y - Y_hat
D1 = matrix(NA,4,4)

for(i in 1:4){
  for(j in 1:4){
    D1[i,j] = sum((R[i,]-R[j,])^2)    
  }
}

D2 = matrix(NA,4,4)
I = diag(rep(1,n))

for(i in 1:4){
  for(j in 1:4){
    D2[i,j] = (I-H)[i,]%*%G%*%(I-H)[,i] + (I-H)[j,]%*%G%*%(I-H)[,j] - 2*(I-H)[i,]%*%G%*%(I-H)[,j]
  }
}



















####### Fake data


fake.data <- data.frame(clus.raw = rep(1:2, each = 5),
                        clus.corrected = rep(1:2, times = 5)) %>% 
  mutate(Temperature = 10 * clus.raw + rnorm(10))
head(fake.data)

fdata <- gather(fake.data, key = "Clustering", value = "Group", clus.raw:clus.corrected)
head(fdata)

## Temperature against clusters for clusters computed on raw distances
ggplot(fake.data, aes(x = clus.raw, y = Temperature)) + geom_point()
## Temperature against clusters for clusters computed on corrected distances
ggplot(fake.data, aes(x = clus.corrected, y = Temperature)) + geom_point()

ggplot(fdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 1)















# ---
#   title: "22_01_2018"
# output: html_document
# ---
#   
#   ```{r setup, include=FALSE}
# knitr::opts_chunk$set(echo = TRUE)
# ```
# 
# ```{r,echo=FALSE, message=FALSE, warning=FALSE}
# rm(list=objects())
# 
# library(vegan)
# 
# ############################################################ Source custom scripts
# 
# source('~/metaord/utils.R')
# 
# ############################################################ Design
# 
# data.wd <- "/home/rnarci/Bureau/CDDrnarci/Donnees/"
# design = read.table(file=file.path(data.wd, "param_bioadvection.csv"),sep="",header=TRUE)
# rownames(design) <- design$Sample
# 
# ############################################################ Import data
# 
# size_fraction = "0.8-5"
# import_data(size_fraction, samples = rownames(design))
# 
# ############################################################ Subset design
# 
# design <- design[metagenomic_sample, ]
# design <- design[,-c(1,6,7)]
# ```
# 
# # Clustering de taille 2
# 
# ```{r,echo=FALSE, message=FALSE, eval=TRUE}
# D_MDS1 = matrix(NA,dim(jaccard_abundance)[1],dim(jaccard_abundance)[1])
# D_MDS2 = matrix(NA,dim(jaccard_abundance)[1],dim(jaccard_abundance)[1])
# 
# mod1 = capscale(as.dist(jaccard_abundance) ~ Temperature + SSD + SI_temperature + Depth + Silicate + SI_nitrates + Temperature:SI_nitrates + Temperature:SSD, data = design)
# mod2 = capscale(as.dist(jaccard_abundance) ~ Temperature, data = design)
# Y_tilde1 = as.matrix(mod1$CCA$u)
# Y_tilde1 = scale(Y_tilde1)
# Y_tilde2 = as.matrix(mod2$CCA$u)
# Y_tilde2 = scale(Y_tilde2)
# 
# parms = c(design$Temperature,design$SSD,design$SI_temperature,design$Depth,design$Silicate,design$SI_nitrates)
# X1 = matrix(parms,nrow=dim(jaccard_abundance)[1],ncol=6)
# X2 = design$Temperature
# H1 = X1%*%solve(t(X1)%*%X1)%*%t(X1)
# H2 = X2%*%solve(t(X2)%*%X2)%*%t(X2)
# 
# Y_hat1 = H1%*%Y_tilde1 
# Y_hat2 = H2%*%Y_tilde2 
# 
# for(i in 1:dim(jaccard_abundance)[1]){
#   for(j in 1:dim(jaccard_abundance)[1]){
#     D_MDS1[i,j] = sqrt(sum((Y_hat1[i,] - Y_hat1[j,])^2)) 
#     D_MDS2[i,j] = sqrt(sum((Y_hat2[i,] - Y_hat2[j,])^2)) 
#   }
# }
# 
# ## Sans MDS
# 
# D_without_MDS1 = matrix(NA,dim(jaccard_abundance)[1],dim(jaccard_abundance)[1])
# D_without_MDS2 = matrix(NA,dim(jaccard_abundance)[1],dim(jaccard_abundance)[1])
# 
# A = -0.5*jaccard_abundance^2
# Id = diag(identity(dim(jaccard_abundance)[1]))
# un = t(t(rep(1,dim(jaccard_abundance)[1])))
# scale = Id-(1/dim(jaccard_abundance)[1])*un%*%t(un)
# 
# G = scale%*%A%*%scale
# 
# for(i in 1:dim(jaccard_abundance)[1]){
#   for(j in 1:dim(jaccard_abundance)[1]){
#     D_without_MDS1[i,j] = sqrt(H1[i,]%*%G%*%H1[,i] + H1[j,]%*%G%*%H1[,j] - 2*H1[i,]%*%G%*%H1[,j])
#     D_without_MDS2[i,j] = sqrt(H2[i,]%*%G%*%H2[,i] + H2[j,]%*%G%*%H2[,j] - 2*H2[i,]%*%G%*%H2[,j])
#   }
# }
# 
# ####### Pre-requis pour clustering par k-means
# 
# fit1 = cmdscale(jaccard_abundance, eig=TRUE, k=2)
# fit2 = cmdscale(D_MDS1, eig=TRUE, k=2)
# fit3 = cmdscale(D_without_MDS1, eig=TRUE, k=2)
# fit4 = cmdscale(D_MDS2, eig=TRUE, k=2)
# fit5 = cmdscale(D_without_MDS2, eig=TRUE, k=2)
# 
# ####### Pre-requis pour clustering ordinal
# 
# library(loe)
# library(fcd)
# library(clues)
# 
# kNN1 = make.kNNG(jaccard_abundance, k = 35, symm = TRUE, weight = FALSE)
# kNN2 = make.kNNG(D_MDS1, k = 35, symm = TRUE, weight = FALSE)
# kNN3 = make.kNNG(D_without_MDS1, k = 35, symm = TRUE, weight = FALSE)
# kNN4 = make.kNNG(D_MDS2, k = 35, symm = TRUE, weight = FALSE)
# kNN5 = make.kNNG(D_without_MDS2, k = 35, symm = TRUE, weight = FALSE)
# 
# ####### Clustering de l clusters
# 
# l = 2
# 
# res1 = kmeans(fit1$points,l)$cluster
# res2 = kmeans(fit2$points,l)$cluster
# res3 = kmeans(fit3$points,l)$cluster
# res4 = kmeans(fit4$points,l)$cluster
# res5 = kmeans(fit5$points,l)$cluster
# 
# res6 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res7 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res8 = spectral.clustering(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res9 = spectral.clustering(kNN4, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res10 = spectral.clustering(kNN5, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# 
# ####### Comparaison des clusterings
# 
# tab = matrix(NA,10,10)
# 
# for(i in 1:10){
#   for(j in 1:10){
#     tab[i,j] = adjustedRand(get(paste("res",i,sep="")),get(paste("res",j,sep="")),randMethod="Rand") 
#   }
# }
# round(tab,2)
# ```
# 
# # Clustering de taille 3
# 
# ```{r,echo=FALSE, message=FALSE, eval=TRUE}
# l = 3
# 
# res1 = kmeans(fit1$points,l)$cluster
# res2 = kmeans(fit2$points,l)$cluster
# res3 = kmeans(fit3$points,l)$cluster
# res4 = kmeans(fit4$points,l)$cluster
# res5 = kmeans(fit5$points,l)$cluster
# 
# res6 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res7 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res8 = spectral.clustering(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res9 = spectral.clustering(kNN4, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res10 = spectral.clustering(kNN5, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# 
# ####### Comparaison des clusterings
# 
# tab = matrix(NA,10,10)
# 
# for(i in 1:10){
#   for(j in 1:10){
#     tab[i,j] = adjustedRand(get(paste("res",i,sep="")),get(paste("res",j,sep="")),randMethod="Rand") 
#   }
# }
# round(tab,2)
# ```
# 
# # Clustering de taille 4
# 
# ```{r,echo=FALSE, message=FALSE, eval=TRUE}
# l = 4
# 
# res1 = kmeans(fit1$points,l)$cluster
# res2 = kmeans(fit2$points,l)$cluster
# res3 = kmeans(fit3$points,l)$cluster
# res4 = kmeans(fit4$points,l)$cluster
# res5 = kmeans(fit5$points,l)$cluster
# 
# res6 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res7 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res8 = spectral.clustering(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res9 = spectral.clustering(kNN4, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# res10 = spectral.clustering(kNN5, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# 
# ####### Comparaison des clusterings
# 
# tab = matrix(NA,10,10)
# 
# for(i in 1:10){
#   for(j in 1:10){
#     tab[i,j] = adjustedRand(get(paste("res",i,sep="")),get(paste("res",j,sep="")),randMethod="Rand") 
#   }
# }
# round(tab,2)
# ```
# 
# ```{r}
# datatable(results) %>% formatRound(columns = 1:7, digits = 3)
# ```
# 
# 
# 
# 
# 
# 
# 
# 
# 




########## Plot geographic map 

# Premiere facon

library(rworldmap)
newmap <- getMap(resolution = "li")
plot(newmap, xlim = c(-180, 90), ylim = c(-75, 75), asp = 1)
points(gps$Mean_longitude, gps$Mean_latitude, col = res4 + 1, cex = .6, pch = 20)

# Deuxieme facon 

sbbox <- make_bbox(lon = gps$lon, lat = gps$lat, f = .1)
sq_map <- get_map(location = sbbox, maptype = "watercolor", source = "google")


ggmap(sq_map) + geom_point(data = gps, mapping = aes(x = gps$lon, y = gps$lat), color = "red") +
  #geom_text(data = gps, aes(label = rownames(gps)), angle = 60, hjust = 0, color = "yellow")

  
  
  
  
  
  
  
  
  



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







































##################################### Quelques tests sur la profondeur : Depth

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

############################################################ Lagrangian distances

lagrangien = read.table(file=file.path(data.wd,"tarrive_min_surface_1000.csv"),header=TRUE, sep="")
colnames(lagrangien) <- rownames(lagrangien)

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

lagrangien = lagrangien[metagenomic_sample,metagenomic_sample]

library(mice)
mice = complete(mice(design,method="norm.predict",m=1))
design$Phosphates = mice$Phosphates
design$NO2NO3 = mice$NO2NO3

design[,c(12,13)] = c(gps$Mean_latitude,gps$Mean_longitude)
colnames(design)[c(12,13)] = c("Latitude","Longitude")

library(loe)
library(fcd)
library(clues)

D = jaccard_abundance
n = dim(D)[1]

parms1 = c(design$Depth)
X1 = matrix(parms1,nrow=n,ncol=1)

H1 = X1%*%solve(t(X1)%*%X1)%*%t(X1)

D_without_MDS1 = matrix(NA,n,n)

J = diag(rep(1,n)) - matrix(1,n,n)/n
G = -0.5*J%*%D^2%*%J

for(i in 1:n){
  for(j in 1:n){
    D_without_MDS1[i,j] = sqrt(H1[i,]%*%G%*%H1[,i] + H1[j,]%*%G%*%H1[,j] - 2*H1[i,]%*%G%*%H1[,j])
  }
}

l = 2

cluster_ref = c()
  
for(i in 1:n){
  if(str_detect(metagenomic_sample,"DCM")[i]==TRUE){
    cluster_ref[i] = 1
  }else
  {
    cluster_ref[i] = 2
  }
}

library(fossil)

# result_cmdscale = matrix(NA,48,2)
# 
# for(k1 in 3:50){
#   fit2 = cmdscale(D_without_MDS1, eig=TRUE, k=k1)
#   res2 = kmeans(fit2$points,l, nstart = 1000)$cluster
#   result_cmdscale[k1-2,1] = rand.index(res2,cluster_ref) #= round(abs(sum(str_detect(metagenomic_sample,"DCM")[res2==1]) - sum(str_detect(metagenomic_sample,"DCM")[res2==2]))/sum(str_detect(metagenomic_sample,"DCM"))*100)
# 
#   fit2bis = cmdscale(D, eig=TRUE, k=k1)
#   res2bis = kmeans(fit2bis$points,l, nstart = 1000)$cluster
#   result_cmdscale[k1-2,2] = rand.index(res2bis,cluster_ref)# = round(abs(sum(str_detect(metagenomic_sample,"DCM")[res2bis==1]) - sum(str_detect(metagenomic_sample,"DCM")[res2bis==2]))/sum(str_detect(metagenomic_sample,"DCM"))*100)
# 
#   #cat(sprintf("Itération %s \n",k1-2))
# }
# round(result_cmdscale*100)
# # mean(result_cmdscale[,1])
# # mean(result_cmdscale[,2])

result_kNNG = matrix(NA,58,2)

for(k2 in 3:60){
  kNN2 = make.kNNG(D_without_MDS1, k = k2, symm = TRUE, weight = FALSE)
  res5 = spectral.clustering.new(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  result_kNNG[k2-2,1] = rand.index(res5,cluster_ref)#= round(abs(sum(str_detect(metagenomic_sample,"DCM")[res5==1]) - sum(str_detect(metagenomic_sample,"DCM")[res5==2]))/sum(str_detect(metagenomic_sample,"DCM"))*100)

  kNN2bis = make.kNNG(D, k = k2, symm = TRUE, weight = FALSE)
  res5bis = spectral.clustering.new(kNN2bis, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  result_kNNG[k2-2,2] = rand.index(res5bis,cluster_ref)#= round(abs(sum(str_detect(metagenomic_sample,"DCM")[res5bis==1]) - sum(str_detect(metagenomic_sample,"DCM")[res5bis==2]))/sum(str_detect(metagenomic_sample,"DCM"))*100)

  cat(sprintf("Itération %s \n",k2-2))
}
round(result_kNNG*100)
# mean(result_kNNG[,1])
# mean(result_kNNG[,2])





############################################### Test de stabilite

rm(list=objects())

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(kernlab)
library(fossil)

############################################################ Source custom scripts

source('~/metaord/utils.R')

############################################################ Design

data.wd <- "/home/rnarci/Bureau/CDDrnarci/Donnees/"
design = read.table(file=file.path(data.wd, "param_bioadvection.csv"),sep="",header=TRUE)
rownames(design) <- design$Sample

############################################################ Import longitude and latitude of the stations

gps = read.table(file=file.path(data.wd,"GPScoordinates2.csv"),header=TRUE,sep="")
rownames(gps) <- gps$Station

############################################################ Lagrangian distances

lagrangien = read.table(file=file.path(data.wd,"tarrive_min_surface_1000.csv"),header=TRUE, sep="")
colnames(lagrangien) <- rownames(lagrangien)

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

lagrangien = lagrangien[metagenomic_sample,metagenomic_sample]

library(mice)
mice = complete(mice(design,method="norm.predict",m=1))
design$Phosphates = mice$Phosphates
design$NO2NO3 = mice$NO2NO3

design[,c(12,13)] = c(gps$Mean_latitude,gps$Mean_longitude)
colnames(design)[c(12,13)] = c("Latitude","Longitude")

library(loe)
library(fcd)
library(clues)

D = jaccard_abundance
n = dim(D)[1]

l = 10

kNN = make.kNNG(D, k = 37, symm = TRUE, weight = FALSE)

res1 = specc(x = kNN, centers = l, iterations = 1000)
res2 = spectral.clustering(kNN, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
rand.index(res1,res2)

res1 = spectral.clustering(kNN, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res2 = spectral.clustering(kNN, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
rand.index(res1,res2)

res1 = specc(x = kNN, centers = l, iterations = 200)
res2 = specc(x = kNN, centers = l, iterations = 200)
rand.index(res1,res2)

fit = cmdscale(D,eig=T,k=n-1)
res1 = kmeans(fit$points, centers = l, nstart = 1000)$cluster
res2 = kmeans(fit$points, centers = l, nstart = 1000)$cluster
rand.index(res1,res2)

res1 = spectral.clustering.new(kNN, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res2 = spectral.clustering.new(kNN, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
rand.index(res1,res2)



#############################" Tests autres matrices de distances et fractions de taille (05/03/2018)

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

############################################################ Lagrangian distances

lagrangien = read.table(file=file.path(data.wd,"tarrive_min_surface_1000.csv"),header=TRUE, sep="")
colnames(lagrangien) <- rownames(lagrangien)

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

lagrangien = lagrangien[metagenomic_sample,metagenomic_sample]

library(mice)
mice = complete(mice(design,method="norm.predict",m=1))
design$Phosphates = mice$Phosphates
design$NO2NO3 = mice$NO2NO3

design[,c(12,13)] = c(gps$Mean_latitude,gps$Mean_longitude)
colnames(design)[c(12,13)] = c("Latitude","Longitude")

D = jaccard_abundance
n = dim(D)[1]




parms1 = design$Depth
X1 = matrix(parms1,nrow=n,ncol=1)


H1 = X1%*%solve(t(X1)%*%X1)%*%t(X1)

D_without_MDS1 = matrix(NA,n,n)


J = diag(rep(1,n)) - matrix(1,n,n)/n
G = -0.5*J%*%D^2%*%J

for(i in 1:n){
  for(j in 1:n){
    D_without_MDS1[i,j] = sqrt(H1[i,]%*%G%*%H1[,i] + H1[j,]%*%G%*%H1[,j] - 2*H1[i,]%*%G%*%H1[,j])
  }
}


fit2 = cmdscale(D_without_MDS1, eig=TRUE, k=n-1)




library(loe)
library(fcd)
library(clues)


kNN2 = make.kNNG(D_without_MDS1, k = 37, symm = TRUE, weight = FALSE)


####### Clustering de l clusters

l = 2


res2 = kmeans(fit2$points,l, nstart = 1000)$cluster



res5 = spectral.clustering.new(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)


library(rworldmap)
newmap <- getMap(resolution = "li")

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


color = c()
labels = c()
cl = res5

for(i in 1:length(cl)){
  if(cl[i]==1){
    color[i] = "blue"
  }
  if(cl[i]==2){
    color[i] = "red"
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



########################## Stocker les clusterings pour chaque matrice, nombre de clusters et nombre 
########################## de plus proches voisins

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


########################################## Representation des resultats pour stability_study

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

p = ggplot(res, aes(x = Diff.k, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)
ggsave(p, file = "all_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/RI_spectral", width = 15, height = 10, dpi = 100)

p = ggplot(res %>% filter(k1 > 10 & k2 > 10), aes(x = Diff.k, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)
ggsave(p, file = "seuil1_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/RI_spectral", width = 15, height = 10, dpi = 100)

p = ggplot(res %>% filter(k1 > 30 & k2 > 30), aes(x = Diff.k, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)
ggsave(p, file = "seuil2_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/RI_spectral", width = 15, height = 10, dpi = 100)

p = ggplot(res %>% filter(Diff.k <= 20), aes(x = Diff.k, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)
ggsave(p, file = "seuil3_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/RI_spectral", width = 15, height = 10, dpi = 100)

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

p = ggplot(res, aes(x = Diff.d, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)
ggsave(p, file = "all_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/RI_kmeans", width = 15, height = 10, dpi = 100)

p = ggplot(res %>% filter(d1 > 10 & d2 > 10), aes(x = Diff.d, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)
ggsave(p, file = "seuil1_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/RI_kmeans", width = 15, height = 10, dpi = 100)

p = ggplot(res %>% filter(d1 > 30 & d2 > 30), aes(x = Diff.d, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)
ggsave(p, file = "seuil2_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/RI_kmeans", width = 15, height = 10, dpi = 100)

p = ggplot(res %>% filter(Diff.d <= 20), aes(x = Diff.d, y = RI, color = Distance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # geom_density_2d() +
  facet_grid(Distance~Nb_Cluster)
ggsave(p, file = "seuil3_points.png",  path="~/Bureau/CDDrnarci/Reunions/14_03_2018/RI_kmeans", width = 15, height = 10, dpi = 100)

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
  scale_y_continuous(limits = c(0,1)) +
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
  scale_y_continuous(limits = c(0,1)) +
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
  scale_y_continuous(limits = c(0,1)) +
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
  scale_y_continuous(limits = c(0,1)) +
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

# oldwd <- getwd() 
# 
# setwd(dir="~/Bureau/CDDrnarci/Donnees")
# data1 = read.xls(xls="Genocenoses_env_parameters_all_tara.xlsx", method = "csv")
# data.0_0.2 = data1[which((as.character(data1$Fraction) == "0-0.2") == TRUE),c(1,3)] 
# data.0.22_3 = data1[which((as.character(data1$Fraction) == "0.22-3") == TRUE),c(1,3)] 
# data.0.8_5 = data1[which((as.character(data1$Fraction) == "0.8-5") == TRUE),c(1,3)] 
# data.5_20 = data1[which((as.character(data1$Fraction) == "May-20") == TRUE),c(1,3)] 
# data.20_180 = data1[which((as.character(data1$Fraction) == "20-180") == TRUE),c(1,3)] 
# data.180_2000 = data1[which((as.character(data1$Fraction) == "180-2000") == TRUE),c(1,3)] 
# 
# rownames(data.0_0.2) = data.0_0.2$Station
# rownames(data.0.22_3) = data.0.22_3$Station
# rownames(data.0.8_5) = data.0.8_5$Station
# rownames(data.5_20) = data.5_20$Station
# rownames(data.20_180) = data.20_180$Station
# rownames(data.180_2000) = data.180_2000$Station
# 
# setwd(oldwd)

for(size_fraction in c("0-0.2","0.22-3","0.8-5","5-20","20-180","180-2000")){
  # if(size_fraction=="0.8-5"){
  #   import_data(size_fraction, samples = rownames(data.0.8_5))
  # }
  # 
  # if(size_fraction=="5-20"){
  #   import_data(size_fraction, samples = rownames(data.5_20))
  # }
  # 
  # if(size_fraction=="180-2000"){
  #   import_data(size_fraction, samples = rownames(data.180_2000))
  # }
  # 
  # if(size_fraction=="0.22-3"){
  #   import_data(size_fraction, samples = rownames(data.0.22_3))
  # }
  # 
  # if(size_fraction=="0-0.2"){
  #   import_data(size_fraction, samples = rownames(data.0_0.2))
  # }
  # 
  # if(size_fraction=="20-180"){
  #   import_data(size_fraction, samples = rownames(data.20_180))
  # }
  # n = dim(jaccard_abundance)[1]
  n = 300
  method = "RI_plot"
  type = "kmeans"
  
  threshold = NULL
  p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
  ggsave(p, file = "RI_plot_all_kmeans.png", width = 15, height = 10, dpi = 100)   
  
  threshold = c(11,n-1)
  p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
  ggsave(p, file = "RI_plot_seuil1_kmeans.png", width = 15, height = 10, dpi = 100)  
  
  threshold = c(31,n-1)
  p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
  ggsave(p, file = "RI_plot_seuil2_kmeans.png", width = 15, height = 10, dpi = 100)  
  
  threshold = 20
  p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
  ggsave(p, file = "RI_plot_seuil3_kmeans.png", width = 15, height = 10, dpi = 100)  
  
  type = "spectral"
  
  threshold = NULL
  p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
  ggsave(p, file = "RI_plot_all_spectral.png", width = 15, height = 10, dpi = 100)   
  
  threshold = c(11,n-1)
  p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
  ggsave(p, file = "RI_plot_seuil1_spectral.png", width = 15, height = 10, dpi = 100)  
  
  threshold = c(31,n-1)
  p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
  ggsave(p, file = "RI_plot_seuil2_spectral.png", width = 15, height = 10, dpi = 100)  
  
  threshold = 20
  p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
  ggsave(p, file = "RI_plot_seuil3_spectral.png", width = 15, height = 10, dpi = 100) 
  
  method = "RI_summary"
  type = NULL
  
  threshold = NULL
  p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
  ggsave(p, file = "RI_summary_all.png", width = 15, height = 10, dpi = 100)   
  
  threshold = c(11,n-1)
  p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
  ggsave(p, file = "RI_summary_seuil1.png", width = 15, height = 10, dpi = 100)  
  
  threshold = c(31,n-1)
  p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
  ggsave(p, file = "RI_summary_seuil2.png", width = 15, height = 10, dpi = 100)  
  
  threshold = 20
  p = stability_study_results(type = type, method = method, threshold = threshold, size_fraction = size_fraction)
  ggsave(p, file = "RI_summary_seuil3.png", width = 15, height = 10, dpi = 100) 
}






################################ Recherche artisanale des differences entre K-means et hclust

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

comparison = compare_clusterings()

size_fraction = "0-0.2"
dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0-0.2/mat_kmeans")
setwd(dir=dir)
l = max(data.0_0.2$Genocenose)
index = which(comparison$Distance=="SBP" & comparison$Fraction==size_fraction & comparison$Type=="Kmeans")
new_RI = comparison[index,]$RI
k = which.max(new_RI)
RI_max = new_RI[k]
mat = read.table("Mat11_kmeans")

h_clust = data.0_0.2$Genocenose
n = dim(data.0_0.2)[1]
K_means = rep(NA,n)

for(i in 3:(n+2)){
  K_means[i-2] = mat[which(mat[,1]==k & mat[,2]==l),i]
}
mat[which(mat[,1]==k & mat[,2]==l),3:(n+2)]-K_means
rand.index(K_means,h_clust)==RI_max

K_means
h_clust
K_means_ex = K_means
h_clust_ex = h_clust
K_means[which(K_means_ex==6)]=7
K_means[which(K_means_ex==7)]=3
K_means[which(K_means_ex==3)]=6
K_means
h_clust
h_clust[which(h_clust_ex==4)]=2
h_clust[which(h_clust_ex==2)]=8
h_clust[which(h_clust_ex==8)]=4
h_clust_ex_ex = h_clust_ex
h_clust_ex = h_clust
h_clust[which(h_clust_ex==4)]=5
h_clust[which(h_clust_ex==5)]=4
K_means
h_clust
rand.index(h_clust_ex_ex,K_means_ex)
rand.index(h_clust,K_means)
rand.index(K_means,K_means_ex)
rand.index(h_clust,h_clust_ex_ex)

rownames(data.0_0.2)[h_clust==2]
rownames(data.0_0.2)[K_means==1]

rownames(data.0_0.2)[K_means==8]
rownames(data.0_0.2)[h_clust==1]

rownames(data.0_0.2)[K_means==4]
rownames(data.0_0.2)[h_clust==6]

#############################################################################################

size_fraction = "0.22-3"
dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.22-3/mat_kmeans")
setwd(dir=dir)
l = max(data.0.22_3$Genocenose)
index = which(comparison$Distance=="SBP" & comparison$Fraction==size_fraction & comparison$Type=="Kmeans")
new_RI = comparison[index,]$RI
k = which.max(new_RI)
RI_max = new_RI[k]
which(matrices.list=="sorensen_braycurtis_prevalence")
mat = read.table("Mat11_kmeans")

h_clust = data.0.22_3$Genocenose
n = dim(data.0.22_3)[1]
K_means = rep(NA,n)

for(i in 3:(n+2)){
  K_means[i-2] = mat[which(mat[,1]==k & mat[,2]==l),i]
}
mat[which(mat[,1]==k & mat[,2]==l),3:(n+2)]-K_means
rand.index(K_means,h_clust)==RI_max

K_means
h_clust
K_means_ex = K_means
h_clust_ex = h_clust
K_means[which(K_means_ex==1)]=5
K_means[which(K_means_ex==5)]=1
K_means
h_clust
h_clust[which(h_clust_ex==7)]=6
h_clust[which(h_clust_ex==6)]=2
h_clust[which(h_clust_ex==2)]=7
h_clust[which(h_clust_ex==3)]=8
h_clust[which(h_clust_ex==8)]=3
K_means
h_clust
rand.index(h_clust_ex,K_means_ex)
rand.index(h_clust,K_means)
rand.index(K_means,K_means_ex)
rand.index(h_clust,h_clust_ex)

rownames(data.0.22_3)[h_clust==6]
rownames(data.0.22_3)[K_means==6]

rownames(data.0.22_3)[K_means==5]
rownames(data.0.22_3)[h_clust==5]

rownames(data.0.22_3)[K_means==8]
rownames(data.0.22_3)[h_clust==8]

#############################################################################################

size_fraction = "0.8-5"
dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.8-5/mat_kmeans")
setwd(dir=dir)
l = max(data.0.8_5$Genocenose)
index = which(comparison$Distance=="SBP" & comparison$Fraction==size_fraction & comparison$Type=="Kmeans")
new_RI = comparison[index,]$RI
k = which.max(new_RI)
RI_max = new_RI[k]
which(matrices.list=="sorensen_braycurtis_prevalence")
mat = read.table("Mat13_kmeans")

h_clust = data.0.8_5$Genocenose
n = dim(data.0.8_5)[1]
K_means = rep(NA,n)

for(i in 3:(n+2)){
  K_means[i-2] = mat[which(mat[,1]==k & mat[,2]==l),i]
}
mat[which(mat[,1]==k & mat[,2]==l),3:(n+2)]-K_means
rand.index(K_means,h_clust)==RI_max

K_means
h_clust
K_means_ex = K_means
h_clust_ex = h_clust
K_means[which(K_means_ex==11)]=9
K_means[which(K_means_ex==9)]=11
K_means[which(K_means_ex==2)]=8
K_means[which(K_means_ex==8)]=2
K_means[which(K_means_ex==1)]=10
K_means[which(K_means_ex==10)]=1
K_means
h_clust
h_clust[which(h_clust_ex==11)]=4
h_clust[which(h_clust_ex==4)]=11
h_clust[which(h_clust_ex==7)]=2
h_clust[which(h_clust_ex==2)]=7
K_means
h_clust
rand.index(h_clust_ex,K_means_ex)
rand.index(h_clust,K_means)
rand.index(K_means,K_means_ex)
rand.index(h_clust,h_clust_ex)

rownames(data.0.8_5)[h_clust==9]
rownames(data.0.8_5)[K_means==9]

rownames(data.0.8_5)[K_means==8]
rownames(data.0.8_5)[h_clust==8]

rownames(data.0.8_5)[K_means==3]
rownames(data.0.8_5)[h_clust==3]

rownames(data.0.8_5)[K_means==4]
rownames(data.0.8_5)[h_clust==4]

rownames(data.0.8_5)[K_means==7]
rownames(data.0.8_5)[h_clust==7]

#############################################################################################

size_fraction = "5-20"
dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/5-20/mat_kmeans")
setwd(dir=dir)
l = max(data.5_20$Genocenose)
index = which(comparison$Distance=="SBP" & comparison$Fraction==size_fraction & comparison$Type=="Kmeans")
new_RI = comparison[index,]$RI
k = which.max(new_RI)
RI_max = new_RI[k]
which(matrices.list=="sorensen_braycurtis_prevalence")
mat = read.table("Mat13_kmeans")

h_clust = data.5_20$Genocenose
n = dim(data.5_20)[1]
K_means = rep(NA,n)

for(i in 3:(n+2)){
  K_means[i-2] = mat[which(mat[,1]==k & mat[,2]==l),i]
}
mat[which(mat[,1]==k & mat[,2]==l),3:(n+2)]-K_means
rand.index(K_means,h_clust)==RI_max

K_means
h_clust
K_means_ex = K_means
h_clust_ex = h_clust

rand.index(h_clust_ex,K_means_ex)
rand.index(h_clust,K_means)
rand.index(K_means,K_means_ex)
rand.index(h_clust,h_clust_ex)

rownames(data.5_20)[h_clust==4]
rownames(data.5_20)[K_means==6]
rownames(data.5_20)[K_means==1]

rownames(data.5_20)[K_means==4]
rownames(data.5_20)[K_means==3]
rownames(data.5_20)[h_clust==3]

rownames(data.5_20)[K_means==5]
rownames(data.5_20)[K_means==6]
rownames(data.5_20)[h_clust==5]
rownames(data.5_20)[h_clust==6]

#############################################################################################

size_fraction = "20-180"
dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/20-180/mat_kmeans")
setwd(dir=dir)
l = max(data.20_180$Genocenose)
index = which(comparison$Distance=="SBP" & comparison$Fraction==size_fraction & comparison$Type=="Kmeans")
new_RI = comparison[index,]$RI
k = which.max(new_RI)
RI_max = new_RI[k]
which(matrices.list=="sorensen_braycurtis_prevalence")
mat = read.table("Mat11_kmeans")

h_clust = data.20_180$Genocenose
n = dim(data.20_180)[1]
K_means = rep(NA,n)

for(i in 3:(n+2)){
  K_means[i-2] = mat[which(mat[,1]==k & mat[,2]==l),i]
}
mat[which(mat[,1]==k & mat[,2]==l),3:(n+2)]-K_means
rand.index(K_means,h_clust)==RI_max

K_means
h_clust
K_means_ex = K_means
h_clust_ex = h_clust
K_means[which(K_means_ex==3)]=6
K_means[which(K_means_ex==6)]=3
K_means
h_clust
h_clust[which(h_clust_ex==4)]=3
h_clust[which(h_clust_ex==3)]=4
K_means
h_clust
rand.index(h_clust_ex,K_means_ex)
rand.index(h_clust,K_means)
rand.index(K_means,K_means_ex)
rand.index(h_clust,h_clust_ex)

rownames(data.20_180)[h_clust==6]
rownames(data.20_180)[K_means==6]

rownames(data.20_180)[K_means==5]
rownames(data.20_180)[K_means==2]
rownames(data.20_180)[h_clust==5]

rownames(data.20_180)[K_means==3]
rownames(data.20_180)[h_clust==3]

#############################################################################################

size_fraction = "180-2000"
dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/180-2000/mat_kmeans")
setwd(dir=dir)
l = max(data.180_2000$Genocenose)
index = which(comparison$Distance=="SBP" & comparison$Fraction==size_fraction & comparison$Type=="Kmeans")
new_RI = comparison[index,]$RI
k = which.max(new_RI)
RI_max = new_RI[k]
which(matrices.list=="sorensen_braycurtis_prevalence")
mat = read.table("Mat13_kmeans")

h_clust = data.180_2000$Genocenose
n = dim(data.180_2000)[1]
K_means = rep(NA,n)

for(i in 3:(n+2)){
  K_means[i-2] = mat[which(mat[,1]==k & mat[,2]==l),i]
}
mat[which(mat[,1]==k & mat[,2]==l),3:(n+2)]-K_means
rand.index(K_means,h_clust)==RI_max

K_means
h_clust
K_means_ex = K_means
h_clust_ex = h_clust
K_means[which(K_means_ex==6)]=8
K_means[which(K_means_ex==8)]=6
K_means
h_clust
rand.index(h_clust_ex,K_means_ex)
rand.index(h_clust,K_means)
rand.index(K_means,K_means_ex)
rand.index(h_clust,h_clust_ex)

rownames(data.0.8_5)[h_clust==8]
rownames(data.0.8_5)[K_means==8]

rownames(data.0.8_5)[K_means==4]
rownames(data.0.8_5)[K_means==3]
rownames(data.0.8_5)[K_means==2]
rownames(data.0.8_5)[h_clust==5]













################################################ Etude des valeurs propres pour le MDS

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

size_fraction = "0-0.2"
delete = which(is.na(data.0_0.2$Genocenose))
data.0_0.2 = data.0_0.2[-delete,]
import_data(size_fraction, samples = rownames(data.0_0.2))
data.0_0.2 = data.0_0.2[metagenomic_sample,]
n = dim(jaccard_abundance)[1]
l = 8

mat = "ochiai_abundance"
seuil = 25

test = list()
points = list()
clust = list()
ri = c()
equal_vp = c()

for(k in 1:round(n/2)){
  mds = cmdscale(get(mat), k = k, eig = T)
  test[[k]] = mds$eig
  points[[k]] = mds$points
}

for(k in 1:round(n/2)){
  equal_vp[k] = sum(test[[seuil]]==test[[k]])
}

for(k in 1:round(n/2)){
  clust[[k]] = kmeans(x = points[[k]], centers = l, nstart = 1000, iter.max = 20)$cluster
}

for(k in 1:round(n/2)){
  ri[k]=rand.index(clust[[seuil]],clust[[k]])
}

equal_vp
ri

######################################### Decroissance des valeurs propres

test[[seuil]]
plot(1:68, test[[seuil]], xlab="Nombre de dimension", ylab="Valeur propre", col = "blue", pch = 19)


























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
library(poilog)
library(vegan)

source('~/metaord/utils.R')

### plot density for given parameters 
barplot(dpoilog(n=0:20,mu=2,sig=1),names.arg=0:20)

### draw random deviates from a community of 50 species 
rpoilog(S=50,mu=2,sig=1)

### draw random deviates including zeros 
rpoilog(S=50,mu=2,sig=1,keep0=TRUE)

### draw random deviates with sampling intensity = 0.5 
rpoilog(S=50,mu=2,sig=1,nu=0.5)

### how many species are likely to be observed 
### (given S,mu,sig2 and nu)? 
hist(replicate(1000,length(rpoilog(S=30,mu=0,sig=3,nu=0.7))))

### how many individuals are likely to be observed
### (given S,mu,sig2 and nu)? 
hist(replicate(1000,sum(rpoilog(S=30,mu=0,sig=3,nu=0.7))))



n = 30
p = 40
count_data = matrix(NA,n,p)
# count_databis = matrix(NA,200,50)

mu1 = runif(p,-1,1)
mu2 = runif(p,-0.9,0.9)
sig1 = 0.5
sig2 = 0.45

# mu3 = runif(50,-1,1)
# 
# sigma1 = rep(1,50)
# sigma2 = rep(1.1,50)
# sigma3 = rep(1,50)
# 
sum(mu1 + sig1^2/2)
sum(mu2 + sig2^2/2)
# sum(mu3 + sigma3^2/2)
# 
# Z1 = rnorm(50, mean = mu1, sd = sigma1)
# Z2 = rnorm(50, mean = mu2, sd = sigma2)
# group1 = t(replicate(100,rpois(50, exp(Z1))))
# group2 = t(replicate(100,rpois(50, exp(Z2))))

group1 = t(replicate(n/2,rpoilog(S=p,mu=mu1,sig=sig1,keep0=TRUE)))
group2 = t(replicate(n/2,rpoilog(S=p,mu=mu2,sig=sig2,keep0=TRUE)))

count_data = rbind(group1,group2)
# rowSums(count_data)
# sum(rowSums(count_data)[1:(n/2)])
# sum(rowSums(count_data)[(n/2+1):n])

braycurtis_distance = as.matrix(vegdist(count_data, method = "bray"))
clust_ref = c(rep(1,n/2),rep(2,n/2))

RI1 = c()
for(k in 1:(n-1)){
  fit1 = cmdscale(braycurtis_distance, k = k)
  clust1 = kmeans(x = fit1, centers = 2, iter.max = 20, nstart = 1000)$cluster
  RI1[k] = rand.index(clust1, clust_ref)
}
RI1

RI2 = c()
for(k in 1:(n-1)){
  kNN1 = make.kNNG(braycurtis_distance, k = k, symm = TRUE, weight = FALSE)
  clust2 = spectral.clustering.new(A = kNN1, K = 2)
  RI2[k] = rand.index(clust2, clust_ref)  
}
RI2

braycurtis_distance = vegdist(count_data, method = "bray")
# adonis(braycurtis_distance ~ clust_ref)

ord = get.order(as.matrix(braycurtis_distance))
embedding = SOE(CM = ord, N = nrow(as.matrix(braycurtis_distance)), p = 7, c = 0.1, maxit = 1000, report = 100)
# summary(lm(embedding$X ~ clust_ref))
# summary(manova(embedding$X ~ clust_ref + runif(n,-1,1)))

RI3 = c()
fit3 = embedding$X
clust3 = kmeans(x = fit3, centers = 2, iter.max = 20, nstart = 1000)$cluster
RI3 = rand.index(clust3, clust_ref)
RI3
