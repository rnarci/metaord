rm(list=objects())

library(kknn)

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

data   = square

par(mfrow=c(1,1))

if(data[1,1]==moons[1,1]){
  
  plot(moons[1:half_n,1],moons[1:half_n,2],xlim=c(-1.5,2),ylim=c(-2.5,1.5),xlab="",ylab="",col="green")
  points(moons[(half_n+1):(2*half_n),1],moons[(half_n+1):(2*half_n),2],col="red")
  
}else{
  
  plot(square[1:half_n,1],square[1:half_n,2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col="green")
  points(square[(half_n+1):(2*half_n),1],square[(half_n+1):(2*half_n),2],col="red")
}

############################## Spectral clustering based on kNN graph directly from data set

l = 2
k = 2

res=specClust(data, centers=l, nn = k)

par(mfrow=c(1,2))

if(data[1,1]==moons[1,1]){
  plot(moons[1:half_n,1],moons[1:half_n,2],xlim=c(-1.5,2),ylim=c(-2.5,1.5),xlab="",ylab="",col="green")
  points(moons[(half_n+1):(2*half_n),1],moons[(half_n+1):(2*half_n),2],col="red")
  
  plot(moons[which(res$cluster==1),1],moons[which(res$cluster==1),2],xlim=c(-1.5,2),ylim=c(-2.5,1.5),xlab="",ylab="",col="red")
  points(moons[which(res$cluster==2),1],moons[which(res$cluster==2),2],col="green")
}else{
  plot(square[1:half_n,1],square[1:half_n,2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col="green")
  points(square[(half_n+1):(2*half_n),1],square[(half_n+1):(2*half_n),2],col="red")
  
  plot(square[which(res$cluster==1),1],square[which(res$cluster==1),2],xlim=c(-1,4),ylim=c(-1,4),xlab="",ylab="",col="red")
  points(square[which(res$cluster==2),1],square[which(res$cluster==2),2],col="green")
}
