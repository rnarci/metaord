rm(list=objects())

library(fcd)

############################## Create data set

half_n = 50

########## 1) Two-moons data set

moons = shapes.two.moon(numObjects=half_n,shape1a=0.4,shape2b=1.2,shape1rFrom=0.8,
                        shape1rTo=1.2,shape2rFrom=0.8, shape2rTo=1.2, outputCsv="", outputCsv2="", 
                        outputColNames=TRUE, outputRowNames=TRUE)
moons = as.matrix(moons$data)

########## 2) Two-squares data set

dim1 = c(runif(half_n,-1,1),runif(half_n,2,4))
dim2 = c(runif(half_n,0,2),runif(half_n,0,2))

square = matrix(c(dim1,dim2),nrow=2*half_n,ncol=2)

########## 3) Choose the data set

data   = moons

par(mfrow=c(1,1))

if(data[1,1]==moons[1,1]){

plot(moons[1:half_n,1],moons[1:half_n,2],xlim=c(-1.5,2),ylim=c(-2.5,1.5),xlab="",ylab="",col="green")
points(moons[(half_n+1):(2*half_n),1],moons[(half_n+1):(2*half_n),2],col="red")

}else{

plot(square[1:half_n,1],square[1:half_n,2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col="green")
points(square[(half_n+1):(2*half_n),1],square[(half_n+1):(2*half_n),2],col="red")
}

############################## Preliminary

########## Compute distance matrix from data set

L = as.matrix(dist(data,method="euclidean",diag=TRUE)) # distance matrix

########## Collect all statements from distance matrix

S = collect_all_statements_distance(L)

########## Optional : collect an arbitrary collection of statements from all statements

n_statements = 10000
u = ceiling(runif(n_statements,0,dim(S)[1]))
S_subset = S[u,]

############################## Perform Algorithm 5 Clustering

n_data = half_n*2
k = 2
l = 2
sigma = 5

algo = k_RNG_clustering_unweighted(S=S_subset,n_data=n_data,k=k,l=l)

kRNG = algo$kRNG 

for(i in 1:(n_data-1)){
  for(j in (i+1):n_data){
    if(kRNG[i,j]==1){
      lines(c(data[i,1],data[j,1]),c(data[i,2],data[j,2]),col="blue")
      }
  }
}

res = algo$res

par(mfrow=c(1,2))

if(data[1,1]==moons[1,1]){
  plot(moons[1:half_n,1],moons[1:half_n,2],xlim=c(-1.5,2),ylim=c(-2.5,1.5),xlab="",ylab="",col="green")
  points(moons[(half_n+1):(2*half_n),1],moons[(half_n+1):(2*half_n),2],col="red")
  
  plot(moons[which(res==1),1],moons[which(res==1),2],xlim=c(-1.5,2),ylim=c(-2.5,1.5),xlab="",ylab="",col="red")
  points(moons[which(res==2),1],moons[which(res==2),2],col="green")
}else{
  plot(square[1:half_n,1],square[1:half_n,2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col="green")
  points(square[(half_n+1):(2*half_n),1],square[(half_n+1):(2*half_n),2],col="red")

  plot(square[which(res==1),1],square[which(res==1),2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col="red")
  points(square[which(res==2),1],square[which(res==2),2],col="green")
}


############################## Clustering with true k-RNG

n_data = half_n*2
k = 2
l = 2

algo = true_kRNG_clustering_unweighted(L=L,k=k,l=l)

kRNG = algo$kRNG 

for(i in 1:(n_data-1)){
  for(j in (i+1):n_data){
    if(kRNG[i,j]==1){
      lines(c(data[i,1],data[j,1]),c(data[i,2],data[j,2]),col="blue")
    }
  }
}

res = algo$res

par(mfrow=c(1,2))

if(data[1,1]==moons[1,1]){
  plot(moons[1:half_n,1],moons[1:half_n,2],xlim=c(-1.5,2),ylim=c(-2.5,1.5),xlab="",ylab="",col="green")
  points(moons[(half_n+1):(2*half_n),1],moons[(half_n+1):(2*half_n),2],col="red")
  
  plot(moons[which(res==1),1],moons[which(res==1),2],xlim=c(-1.5,2),ylim=c(-2.5,1.5),xlab="",ylab="",col="red")
  points(moons[which(res==2),1],moons[which(res==2),2],col="green")
}else{
  plot(square[1:half_n,1],square[1:half_n,2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col="green")
  points(square[(half_n+1):(2*half_n),1],square[(half_n+1):(2*half_n),2],col="red")
  
  plot(square[which(res==1),1],square[which(res==1),2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col="red")
  points(square[which(res==2),1],square[which(res==2),2],col="green")
}

