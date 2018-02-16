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

size_fraction1 = "0-0.2"
size_fraction2 = "0.22-3"
size_fraction3 = "0.8-5"
size_fraction4 = "5-20"
size_fraction5 = "20-180"
size_fraction6 = "180-2000"

import_data(size_fraction3, samples = NULL)

############################################################ First part : spectral clustering on k-NNG

############################## Calculate the best "k" in order to maximise similarities between graphs

res = c()

for(k in 1:50){
  
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
  
  # kNN1 = make.kNNG(jaccard_abundance, k = k, symm = TRUE, weight = TRUE)
  # 
  # kNN2 = make.kNNG(ab_jaccard_abundance, k = k, symm = TRUE, weight = TRUE)
  # 
  # kNN3 = make.kNNG(braycurtis_abundance, k = k, symm = TRUE, weight = TRUE)
  # 
  # kNN4 = make.kNNG(ab_ochiai_abundance, k = k, symm = TRUE, weight = TRUE)
  # 
  # kNN5 = make.kNNG(ab_sorensen_abundance, k = k, symm = TRUE, weight = TRUE)
  # 
  # kNN6 = make.kNNG(simka_jaccard_abundance, k = k, symm = TRUE, weight = TRUE)
  # 
  # kNN7 = make.kNNG(chord_prevalence, k = k, symm = TRUE, weight = TRUE)
  # 
  # kNN8 = make.kNNG(jaccard_prevalence, k = k, symm = TRUE, weight = TRUE)
  # 
  # kNN9 = make.kNNG(kulczynski_prevalence, k = k, symm = TRUE, weight = TRUE)
  # 
  # kNN10 = make.kNNG(ochiai_prevalence, k = k, symm = TRUE, weight = TRUE)
  # 
  # kNN11 = make.kNNG(whittaker_prevalence, k = k, symm = TRUE, weight = TRUE)
  # 
  # kNN12 = make.kNNG(simka_jaccard_prevalence, k = k, symm = TRUE, weight = TRUE)
  
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
# write.table(x=res,file="program1")



restot1=matrix(0,nrow=14,ncol=2)
restot1bis=matrix(0,nrow=14,ncol=2)



############################## Spectral clustering on k-NNG : similarities between clusterings are calculated from ARI.


k = 35


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

for(l in 2:15)
{
  
  ARI1=rep(0,100)
  RI1 = rep(0,100)
  
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
    
    res7 = spectral.clustering.rw(kNN7, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
    
    res8 = spectral.clustering.rw(kNN8, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
    
    res9 = spectral.clustering.rw(kNN9, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
    
    res10 = spectral.clustering.rw(kNN10, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
    
    res11 = spectral.clustering.rw(kNN11, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
    
    res12 = spectral.clustering.rw(kNN12, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
    
    # res1 = specc(kNN7, centers=l)[1:140]
    # 
    # res2 = specc(kNN8,centers=l)[1:140]
    # 
    # res3 = specc(kNN9, centers=l)[1:140]
    # 
    # res4 = specc(kNN10, centers=l)[1:140]
    # 
    # res5 = specc(kNN11, centers=l)[1:140]
    # 
    # res6 = specc(kNN12, centers=l)[1:140]
    
    # res7 = specc(kNN7, centers=l)[1:140]
    # 
    # res8 = specc(kNN8,centers=l)[1:140]
    # 
    # res9 = specc(kNN9, centers=l)[1:140]
    # 
    # res10 = specc(kNN10, centers=l)[1:140]
    # 
    # res11 = specc(kNN11, centers=l)[1:140]
    # 
    # res12 = specc(kNN12, centers=l)[1:140]
    
    
    for(j in 7:11){
      for(m in (j+1):12){
        ARI1[nb_run] = ARI1[nb_run] + adjustedRand(get(paste("res",j,sep="")),get(paste("res",m,sep="")),randMethod="Rand")
        RI1[nb_run] = RI1[nb_run] + rand.index(get(paste("res",j,sep="")),get(paste("res",m,sep="")))
      }
    }
    
    
    ARI1[nb_run] = ARI1[nb_run]/15
    RI1[nb_run] = RI1[nb_run]/15
  }
  
  restot1[l-1,1] = max(ARI1)
  restot1[l-1,2] = max(RI1)
  
}



############################## MDS + k-means


for(l in 2:15){
  ARI1bis=matrix(0,100,10)
  RI1bis = matrix(0,100,10)
  for(k in 1:10){
    
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
    
    
    
    for(nb_run in 1:100)
    {
      # res1 = kmeans(fit1$points,l,nstart = 1000)$cluster
      # 
      # res2 = kmeans(fit2$points,l,nstart = 1000)$cluster
      # 
      # res3 = kmeans(fit3$points,l,nstart = 1000)$cluster
      # 
      # res4 = kmeans(fit4$points,l,nstart = 1000)$cluster
      # 
      # res5 = kmeans(fit5$points,l,nstart = 1000)$cluster
      # 
      # res6 = kmeans(fit6$points,l,nstart = 1000)$cluster
      
      res7 = kmeans(fit7$points,l,nstart = 1000)$cluster
      
      res8 = kmeans(fit8$points,l,nstart = 1000)$cluster
      
      res9 = kmeans(fit9$points,l,nstart = 1000)$cluster
      
      res10 = kmeans(fit10$points,l,nstart = 1000)$cluster
      
      res11 = kmeans(fit11$points,l,nstart = 1000)$cluster
      
      res12 = kmeans(fit12$points,l,nstart = 1000)$cluster
      
      
      
      for(j in 7:11){
        for(m in (j+1):12){
          ARI1bis[nb_run,k] = ARI1bis[nb_run,k] + adjustedRand(get(paste("res",j,sep="")),get(paste("res",m,sep="")),randMethod="Rand")
          RI1bis[nb_run,k] = RI1bis[nb_run,k] + rand.index(get(paste("res",j,sep="")),get(paste("res",m,sep="")))
        }
      }
      
      
      ARI1bis[nb_run,k] = ARI1bis[nb_run,k]/15
      RI1bis[nb_run,k] = RI1bis[nb_run,k]/15
    }
    
    
  }
  restot1bis[l-1,1] = max(ARI1bis)
  restot1bis[l-1,2] = max(RI1bis)
}




############################################################ Second part : density estimate for k-NNG


library(igraph)

############################## Square data

# Graph-based estimator of the number of clusters

half_n = 30
k = 4
pr = 1
# D = 2

# Toy data

n_data = 2*half_n
set.seed(17)
dim1 = c(runif(half_n,-1,1),runif(half_n,2,4))
dim2 = c(runif(half_n,0,2),runif(half_n,0,2))

square = matrix(c(dim1,dim2),nrow=2*half_n,ncol=2)

plot(square[1:half_n,1],square[1:half_n,2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col="green")
points(square[(half_n+1):(2*half_n),1],square[(half_n+1):(2*half_n),2],col="red")

data = square

L = as.matrix(dist(data,method="euclidean",diag=TRUE))

kNN = make.kNNG(L, k = k, symm = FALSE, weight = FALSE)

par(mfrow=c(2,2))

plot(square[1:half_n,1],square[1:half_n,2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col="green")
points(square[(half_n+1):(2*half_n),1],square[(half_n+1):(2*half_n),2],col="red")

for(i in 1:n_data){
  for(j in 1:n_data){
    if(kNN[i,j]==1){
      lines(c(data[i,1],data[j,1]),c(data[i,2],data[j,2]),col="blue")
    }
  }
}

count_neighbors = 0
# c_D = (pi^(D/2))/(factorial(D/2)) # volume of the unit sphere in D dimensions
# a = 1
sub_kNN = kNN
count = n_data
t = -1
# sub_kNNbis = graph_from_adjacency_matrix(kNN, mode = "directed")
# vertices_to_delete = c()

while(count>pr*n_data & t!=n_data){
  count = n_data
  t = t + 1
  for(i in 1:n_data){
    count_neighbors = sum(kNN[i,]) + sum(kNN[,i])
    if((count_neighbors-k)<=t){
      count = count-1
    }
  }
}


for(i in 1:n_data){
  count_neighbors = sum(kNN[i,]) + sum(kNN[,i])
  # k_nearest_neighbor_position = (order(L[i,]))[k+1] # position of the k-nearest neighbor from the vertex i
  # R_k = L[i,k_nearest_neighbor_position] # distance of the k-nearest neighbor from the vertex i
  # density_estimate = (k-1)/(n_data*c_D*R_k^D) # k-NNG density estimate at vertex i
  if((count_neighbors-k)<=t){
    sub_kNN[i,] = rep(0,n_data)
    sub_kNN[,i] = rep(0,n_data)
    # vertices_to_delete[a] = i
    # a = a + 1
  }
}

# sub_kNNbis = delete_vertices(sub_kNNbis,vertices_to_delete)

plot(square[1:half_n,1],square[1:half_n,2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col="green")
points(square[(half_n+1):(2*half_n),1],square[(half_n+1):(2*half_n),2],col="red")

for(i in 1:n_data){
  for(j in 1:n_data){
    if(sub_kNN[i,j]==1){
      lines(c(data[i,1],data[j,1]),c(data[i,2],data[j,2]),col="blue")
    }
  }
}

sub_kNN = graph_from_adjacency_matrix(sub_kNN, mode = "directed")

max_nb_cluster = n_data
max_length_cluster = n_data

res = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
vertices  = c()

t_cluster = na.omit(dfs(sub_kNN,root=1, neimode="all",unreachable=FALSE)$order)
if(length(t_cluster)>1){
  res[1,1:length(t_cluster)] = t_cluster
  vertices = t_cluster 
}

a = 2

for(i in 2:n_data){
  if(length(which((i!=vertices)==TRUE)) == length(vertices)){
    t_cluster = na.omit(dfs(sub_kNN,root=i, neimode="all",unreachable = FALSE)$order)
    if(length(t_cluster)>1){
      res[a,1:length(t_cluster)] = t_cluster
      vertices = c(vertices,t_cluster) 
      a = a + 1
    }
  } 
}


plot(square[1:half_n,1],square[1:half_n,2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col="green")
points(square[(half_n+1):(2*half_n),1],square[(half_n+1):(2*half_n),2],col="red")

for(i in 1:max_nb_cluster){
  test = c()
  test = na.omit(res[i,])
  if(length(test)>1){
    for(s in 1:length(test)){
      for(t in 1:length(test)){
        if(s!=t & kNN[test[s],test[t]]==1){
          lines(c(data[test[s],1],data[test[t],1]),c(data[test[s],2],data[test[t],2]),col="blue")
        }
      }
    }
  }
}

test = c()
color = 1
plot(0,0,xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col=color,type="n")


for(i in 1:max_nb_cluster){
  test = c()
  test = na.omit(res[i,])
  if(length(test)>1){
    points(square[test,1],square[test,2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col=color) 
    color = color + 1
  }
}

res


############################## Gaussian data

# Graph-based estimator of the number of clusters

n_data = 200
k = 40
# pr = 0.9
# D = 2

# Toy data
set.seed(70)
normal = matrix(NA,nrow=n_data,ncol=2)
normal[,1] = rnorm(n_data,0,3)
normal[,2] = rnorm(n_data,0,3)

par(mfrow=c(1,1))
plot(normal[,1],normal[,2],xlab="",ylab="",col="red")


data = normal

L = as.matrix(dist(data,method="euclidean",diag=TRUE))

kNN = make.kNNG(L, k = k, symm = FALSE, weight = FALSE)

par(mfrow=c(2,2))

plot(normal[,1],normal[,2],xlab="",ylab="",col="red")

for(i in 1:n_data){
  for(j in 1:n_data){
    if(kNN[i,j]==1){
      lines(c(data[i,1],data[j,1]),c(data[i,2],data[j,2]),col="blue")
    }
  }
}

count_neighbors = 0
# c_D = (pi^(D/2))/(factorial(D/2)) # volume of the unit sphere in D dimensions
# a = 1
sub_kNN = kNN
count = n_data
t = -1
# sub_kNNbis = graph_from_adjacency_matrix(kNN, mode = "directed")
# vertices_to_delete = c()

while(count>pr*n_data & t!=n_data){
  count = n_data
  t = t + 1
  for(i in 1:n_data){
    count_neighbors = sum(kNN[i,]) + sum(kNN[,i])
    if((count_neighbors-k)<=t){
      count = count-1
    }
  }
}


for(i in 1:n_data){
  count_neighbors = sum(kNN[i,]) + sum(kNN[,i])
  # k_nearest_neighbor_position = (order(L[i,]))[k+1] # position of the k-nearest neighbor from the vertex i
  # R_k = L[i,k_nearest_neighbor_position] # distance of the k-nearest neighbor from the vertex i
  # density_estimate = (k-1)/(n_data*c_D*R_k^D) # k-NNG density estimate at vertex i
  if((count_neighbors-k)<=t){
    sub_kNN[i,] = rep(0,n_data)
    sub_kNN[,i] = rep(0,n_data)
    # vertices_to_delete[a] = i
    # a = a + 1
  }
}

plot(normal[,1],normal[,2],xlab="",ylab="",col="red")

for(i in 1:n_data){
  for(j in 1:n_data){
    if(sub_kNN[i,j]==1){
      lines(c(data[i,1],data[j,1]),c(data[i,2],data[j,2]),col="blue")
    }
  }
}

sub_kNN = graph_from_adjacency_matrix(sub_kNN, mode = "directed")

max_nb_cluster = n_data
max_length_cluster = n_data

res = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
vertices  = c()

t_cluster = na.omit(dfs(sub_kNN,root=1, neimode="all",unreachable=FALSE)$order)
if(length(t_cluster)>1){
  res[1,1:length(t_cluster)] = t_cluster
  vertices = t_cluster 
}

a = 2

for(i in 2:n_data){
  if(length(which((i!=vertices)==TRUE)) == length(vertices)){
    t_cluster = na.omit(dfs(sub_kNN,root=i, neimode="all",unreachable = FALSE)$order)
    if(length(t_cluster)>1){
      res[a,1:length(t_cluster)] = t_cluster
      vertices = c(vertices,t_cluster) 
      a = a + 1
    }
  } 
}


plot(normal[,1],normal[,2],xlab="",ylab="",col="red")

for(i in 1:max_nb_cluster){
  test = c()
  test = na.omit(res[i,])
  if(length(test)>1){
    for(s in 1:length(test)){
      for(t in 1:length(test)){
        if(s!=t & kNN[test[s],test[t]]==1){
          lines(c(data[test[s],1],data[test[t],1]),c(data[test[s],2],data[test[t],2]),col="blue")
        }
      }
    }
  }
}

test = c()
color = 1
plot(normal[,1],normal[,2],xlab="",ylab="",col=color,type="n")


for(i in 1:max_nb_cluster){
  test = c()
  test = na.omit(res[i,])
  if(length(test)>1){
    points(normal[test,1],normal[test,2],xlab="",ylab="",col=color) 
    color = color + 1
  }
}

res




############################## Real data

data1 = jaccard_abundance

data2 = ab_jaccard_abundance

data3 = braycurtis_abundance

data4 = ab_ochiai_abundance

data5 = ab_sorensen_abundance

data6 = simka_jaccard_abundance

data7 = chord_prevalence

data8 = jaccard_prevalence

data9 = kulczynski_prevalence

data10 = ochiai_prevalence

data11 = whittaker_prevalence

data12 = simka_jaccard_prevalence

data = list(data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12)

############################################ Create k-NNG

k = 35

kNN1 = make.kNNG(jaccard_abundance, k = k, symm = FALSE, weight = FALSE)

kNN2 = make.kNNG(ab_jaccard_abundance, k = k, symm = FALSE, weight = FALSE)

kNN3 = make.kNNG(braycurtis_abundance, k = k, symm = FALSE, weight = FALSE)

kNN4 = make.kNNG(ab_ochiai_abundance, k = k, symm = FALSE, weight = FALSE)

kNN5 = make.kNNG(ab_sorensen_abundance, k = k, symm = FALSE, weight = FALSE)

kNN6 = make.kNNG(simka_jaccard_abundance, k = k, symm = FALSE, weight = FALSE)

kNN7 = make.kNNG(chord_prevalence, k = k, symm = FALSE, weight = FALSE)

kNN8 = make.kNNG(jaccard_prevalence, k = k, symm = FALSE, weight = FALSE)

kNN9 = make.kNNG(kulczynski_prevalence, k = k, symm = FALSE, weight = FALSE)

kNN10 = make.kNNG(ochiai_prevalence, k = k, symm = FALSE, weight = FALSE)

kNN11 = make.kNNG(whittaker_prevalence, k = k, symm = FALSE, weight = FALSE)

kNN12 = make.kNNG(simka_jaccard_prevalence, k = k, symm = FALSE, weight = FALSE)

kNN = list(kNN1,kNN2,kNN3,kNN4,kNN5,kNN6,kNN7,kNN8,kNN9,kNN10,kNN11,kNN12)

sub_kNN1 = kNN1

sub_kNN2 = kNN2

sub_kNN3 = kNN3

sub_kNN4 = kNN4

sub_kNN5 = kNN5

sub_kNN6 = kNN6

sub_kNN7 = kNN7

sub_kNN8 = kNN8

sub_kNN9 = kNN9

sub_kNN10 = kNN10

sub_kNN11 = kNN11

sub_kNN12 = kNN12

#sub_kNN = list(sub_kNN1,sub_kNN2,sub_kNN3,sub_kNN4,sub_kNN5,sub_kNN6,sub_kNN7,sub_kNN8,sub_kNN9,sub_kNN10,sub_kNN11,sub_kNN12)


n_data = dim(jaccard_abundance)[1]
range_pr = 0.1 # seq(0.25,0.3,0.01)
# range_t = seq(0,0.2,0.001)
# D = 5

res_tot1 = matrix(NA,nrow=length(range_pr),ncol=2)

res_tot2 = matrix(NA,nrow=length(range_pr),ncol=2)

res_tot3 = matrix(NA,nrow=length(range_pr),ncol=2)

res_tot4 = matrix(NA,nrow=length(range_pr),ncol=2)

res_tot5 = matrix(NA,nrow=length(range_pr),ncol=2)

res_tot6 = matrix(NA,nrow=length(range_pr),ncol=2)

res_tot7 = matrix(NA,nrow=length(range_pr),ncol=2)

res_tot8 = matrix(NA,nrow=length(range_pr),ncol=2)

res_tot9 = matrix(NA,nrow=length(range_pr),ncol=2)

res_tot10 = matrix(NA,nrow=length(range_pr),ncol=2)

res_tot11 = matrix(NA,nrow=length(range_pr),ncol=2)

res_tot12 = matrix(NA,nrow=length(range_pr),ncol=2)

res_tot = list(res_tot1,res_tot2,res_tot3,res_tot4,res_tot5,res_tot6,res_tot7,res_tot8,res_tot9,res_tot10,res_tot11,res_tot12)

pr = 0
iter = 1
# c_D = (pi^(D/2))/(factorial(D/2)) # volume of the unit sphere in D dimensions
# essai = c()

for(pr in range_pr){
  
  sub_kNN = list(sub_kNN1,sub_kNN2,sub_kNN3,sub_kNN4,sub_kNN5,sub_kNN6,sub_kNN7,sub_kNN8,sub_kNN9,sub_kNN10,sub_kNN11,sub_kNN12)
  
  for(num_matrix in 1:12){
    t = -1
    count_neighbors = 0
    count = n_data
    while(count>pr*n_data & t!=n_data){
      count = n_data
      t = t + 1
      for (i in 1:n_data){
        count_neighbors = sum(kNN[[num_matrix]][i,]) + sum(kNN[[num_matrix]][,i])
        # k_nearest_neighbor_position = (order(data[[num_matrix]][i,]))[k+1] # position of the k-nearest neighbor from the vertex i
        # R_k = data[[num_matrix]][i,k_nearest_neighbor_position] # distance of the k-nearest neighbor from the vertex i
        # density_estimate = (k-1)/(n_data*c_D*R_k^D) # k-NNG density estimate at vertex i
        # essai[i] = density_estimate
        if((count_neighbors-k)<=t){
          count = count - 1
        }
      }
    } 
    for(i in 1:n_data){
      count_neighbors = sum(kNN[[num_matrix]][i,]) + sum(kNN[[num_matrix]][,i])
      if((count_neighbors-k)<=t){
        sub_kNN[[num_matrix]][i,] = rep(0,n_data)
        sub_kNN[[num_matrix]][,i] = rep(0,n_data)
      }
    } 
  }
  
  for(num_matrix in 1:12){
    sub_kNN[[num_matrix]] = graph_from_adjacency_matrix(sub_kNN[[num_matrix]], mode = "directed")
  }
  
  max_nb_cluster = n_data
  max_length_cluster = n_data
  
  res1 = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
  
  res2 = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
  
  res3 = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
  
  res4 = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
  
  res5 = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
  
  res6 = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
  
  res7 = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
  
  res8 = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
  
  res9 = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
  
  res10 = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
  
  res11 = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
  
  res12 = matrix(NA,nrow=max_nb_cluster,ncol=max_length_cluster)
  
  res = list(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12)
  
  for(num_matrix in 1:12){
    vertices  = c()
    
    t_cluster = na.omit(dfs(sub_kNN[[num_matrix]],root=1, neimode="all",unreachable=FALSE)$order)
    if(length(t_cluster)>1){
      res[[num_matrix]][1,1:length(t_cluster)] = t_cluster
      vertices = t_cluster 
    }
    
    a = 2
    
    for(i in 2:n_data){
      if(length(which((i!=vertices)==TRUE)) == length(vertices)){
        t_cluster = na.omit(dfs(sub_kNN[[num_matrix]],root=i, neimode="all",unreachable = FALSE)$order)
        if(length(t_cluster)>1){
          res[[num_matrix]][a,1:length(t_cluster)] = t_cluster
          vertices = c(vertices,t_cluster) 
          a = a + 1
        }
      }  
    }
    
  }
  
  for(num_matrix in 1:12){
    if(length(which(is.na(res[[num_matrix]][,1])==FALSE))!=0){
      res_tot[[num_matrix]][iter,] = c(length(which(is.na(res[[num_matrix]][,1])==FALSE)),pr) 
    }
  }
  iter = iter + 1
}

par(mfrow=c(4,3))

for(num_matrix in 1:12){
  plot(0,0,xlim=c(min(range_pr),max(range_pr)),ylim=c(0,10),type="n")
  points(res_tot[[num_matrix]][,2],res_tot[[num_matrix]][,1])
}

