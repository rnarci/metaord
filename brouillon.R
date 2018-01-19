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





































rm(list=objects())

Y = matrix(c(1,2,4,7,8,9,3,4,2,5,6,4,1,5,8,9,5,7,1,2),4,5)
Y = scale(Y)

X = matrix(c(4,5,3,1,7,8,4,9,1,2,5,6),4,3)
H = X%*%solve(t(X)%*%X)%*%t(X)

# Essai 1

v1 = sum(((H%*%Y)[1,])^2)

D = as.matrix(dist(Y,method="euclidean",diag=TRUE))
A = -0.5*D^2
Id = diag(identity(4))
un = t(t(c(1,1,1,1)))
scale = Id-(1/4)*un%*%t(un)
G = scale%*%A%*%scale
v2 = H[1,]%*%G%*%H[,1]

# Essai 2 

Y_hat = H%*%Y 
D1 = matrix(NA,4,4)

for(i in 1:4){
  for(j in 1:4){
    D1[i,j] = sum((Y_hat[i,]-Y_hat[j,])^2)    
  }
}

D2 = matrix(NA,4,4)

for(i in 1:4){
  for(j in 1:4){
    D2[i,j] = H[i,]%*%G%*%H[,i] + H[j,]%*%G%*%H[,j] - 2*H[i,]%*%G%*%H[,j]
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
