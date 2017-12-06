############################## Collect all statements from distance matrix (symetric)

collect_all_statements_distance <- function(L) {
  n = dim(L)[1]
  S = matrix(NA,nrow=n^3,ncol=4) # statements
  statement = 1
  
  for(i in 1:(n-2)){
      for (j in (i+1):(n-1)){
        for(k in (j+1):n){
          if(L[i,j]<L[j,k] & L[i,k]<L[j,k]){
            S[statement,1] = i
            S[statement,2:4] = sort(c(i,j,k))
            statement = statement + 1           
          }
          if(L[j,i]<L[i,k] & L[j,k]<L[i,k]){
            S[statement,1] = j
            S[statement,2:4] = sort(c(i,j,k))
            statement = statement + 1           
          }
          if(L[k,i]<L[i,j] & L[k,j]<L[i,j]){
            S[statement,1] = k
            S[statement,2:4] = sort(c(i,j,k))
            statement = statement + 1           
          }
        }
      }
    }
  return(na.omit(S))
  }


############################## Create two-moons data set

.numObjects<-function(numObjects,numClusters){
  if(length(numObjects)<numClusters){
    numObjects<-rep(numObjects,numClusters)
  }
  if(length(numObjects)>numClusters){
    numObjects<-numObjects[1:numClusters]
  }
  numObjects
}

.toCsv<-function(file,data,klasy,outputRowNames,outputColNames,csv2){
  if(paste(file,"",sep="")!=""){
    if(csv2){
      write.table(cbind(1:dim(data)[1],klasy,data),file=file,sep=";",dec=",",row.names=outputRowNames,col.names=outputColNames)
    }
    else{
      write.table(cbind(1:dim(data)[1],klasy,data),file=file,row.names=outputRowNames,col.names=outputColNames)
    }
  }
}

shapes.two.moon<-function(numObjects,shape1a,shape2b,shape1rFrom, shape1rTo,shape2rFrom, shape2rTo, outputCsv, outputCsv2, outputColNames, outputRowNames){
  lo<-.numObjects(numObjects,2)
  x <- matrix(0, nrow=sum(lo), ncol=2)
  
  for(i in 1:sum(lo)){
    alpha<-runif(1,0,2*pi)
    if(i>lo[1]){
      r=runif(1,shape2rFrom,shape2rTo)
    }
    else{
      r=runif(1,shape1rFrom,shape1rTo)
    }
    x[i,1]<-r*cos(alpha)
    x[i,2]<-r*sin(alpha)
    if(i<=lo[1]){
      x[i,1]=shape1a+abs(x[i,1])
    }
    else{
      x[i,1]=-abs(x[i,1])
      x[i,2]=x[i,2]-shape2b
    }
  }
  data<-x
  klasy<-c(rep(1,lo[1]),rep(2,lo[2]))
  .toCsv(outputCsv,data,klasy,outputColNames, outputRowNames,FALSE)
  .toCsv(outputCsv2,data,klasy, outputColNames, outputRowNames,TRUE)
  list(data=data,clusters=klasy)
}

############################## Algorithm 5 Clustering : unweighted version

k_RNG_clustering_unweighted <- function(S,k,l,n_data){
  n_statements = dim(S)[1]
  n_data = n_data
  N = matrix(0,nrow=n_data,ncol=n_data)
  D = matrix(0,nrow=n_data,ncol=n_data)
  V = matrix(NA,nrow=n_data,ncol=n_data)
  W = matrix(NA,nrow=n_data,ncol=n_data)
  
  for(i in 1:n_data){
    for(j in 1:n_data){
      for(y in 1:n_statements){
        if(i!=j & (S[y,2]==i | S[y,3]==i | S[y,4]==i) & (S[y,2]==j | S[y,3]==j | S[y,4]==j)){
          D[i,j] = D[i,j]+1
          if(S[y,1]!=i & S[y,1]!=j){
            N[i,j]=N[i,j]+1
          }
        }
      }
      if(D[i,j]==0){
        V[i,j] = 1000000 
      }else{
        V[i,j] = N[i,j]/D[i,j]
      }
      if(V[i,j]<k/(n_data-2)){
        W[i,j] = 1 
      }else{
        W[i,j] = 0
      }
    }
  }
  res = spectral.clustering(W, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  return(list(res=res,kRNG=W))
}


############################## Algorithm 5 Clustering : weighted version

k_RNG_clustering_weighted <- function(S,k,l,n_data,sigma){
  n_statements = dim(S)[1]
  n_data = n_data
  sigma = sigma
  N = matrix(0,nrow=n_data,ncol=n_data)
  D = matrix(0,nrow=n_data,ncol=n_data)
  V = matrix(NA,nrow=n_data,ncol=n_data)
  W = matrix(NA,nrow=n_data,ncol=n_data)
  
  for(i in 1:n_data){
    for(j in 1:n_data){
      for(y in 1:n_statements){
        if(i!=j & (S[y,2]==i | S[y,3]==i | S[y,4]==i) & (S[y,2]==j | S[y,3]==j | S[y,4]==j)){
          D[i,j] = D[i,j]+1
          if(S[y,1]!=i & S[y,1]!=j){
            N[i,j]=N[i,j]+1
          }
        }
      }
      if(D[i,j]==0){
        V[i,j] = 1000000 
      }else{
        V[i,j] = N[i,j]/D[i,j]
      }
      if(V[i,j]<k/(n_data-2)){
        W[i,j] = exp(-(V[i,j]^2)/(sigma^2)) 
      }else{
        W[i,j] = 0
      }
    }
  }
  res = spectral.clustering(W, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  return(list(res=res,kRNG=W))
}


############################## Clustering with true k-RNG

# true_kRNG_clustering <- function(L,k,l) {
#   n = dim(L)[1]
#   count = 0
#   kRNG = matrix(NA,nrow=n,ncol=n)
# 
#   for(i in 1:n){
#     for (j in 1:n){
#       for(m in 1:n){
#         if(i!=j & i!=m & j!=m & L[i,j]>max(L[i,m],L[j,m])){
#           count = count+1
#         }
#       }
#       if(count<k & i!=j){
#         kRNG[i,j] = 1
#         count = 0
#       }
#       if(count>=k & i!=j){
#         kRNG[i,j] = 0
#         count = 0
#       }
#       if(i==j){
#         kRNG[i,j] = 0
#       }
#     }
#   }
#   res = spectral.clustering(kRNG, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   return(list(res=res,kRNG=kRNG))
# }

true_kRNG_clustering_unweighted <- function(L,k,l) {
  n = dim(L)[1]
  count = 0
  kRNG = matrix(NA,nrow=n,ncol=n)
  
  for(i in 1:n){
    for (j in 1:n){
      if(i!=j & length(which((L[i,j]>pmax(L[i,-c(i,j)],L[j,-c(i,j)]))==TRUE))<k){
        kRNG[i,j] = 1
      }
      else{
        kRNG[i,j] = 0
      }
    }
  }
  res = spectral.clustering(kRNG, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  return(list(res=res,kRNG=kRNG))
}

true_kRNG_clustering_weighted <- function(L,k,l,sigma) {
  n = dim(L)[1]
  count = 0
  kRNG = matrix(NA,nrow=n,ncol=n)
  
  for(i in 1:n){
    for (j in 1:n){
      N = length(which((L[i,j]>pmax(L[i,-c(i,j)],L[j,-c(i,j)]))==TRUE))
      if(i!=j & N<k){
        kRNG[i,j] = exp(-N^2/((n-2)^2*sigma^2))
      }
      else{
        kRNG[i,j] = 0
      }
    }
  }
  res = spectral.clustering(kRNG, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
  return(list(res=res,kRNG=kRNG))
}

make.kRNG <- function(L,k,sigma){
  n = dim(L)[1]
  count = 0
  kRNG = matrix(NA,nrow=n,ncol=n)
  
  for(i in 1:n){
    for (j in 1:n){
      N = length(which((L[i,j]>pmax(L[i,-c(i,j)],L[j,-c(i,j)]))==TRUE))
      if(i!=j & N<k){
        kRNG[i,j] = 1 #exp(-N^2/((n-2)^2*sigma^2))
      }
      else{
        kRNG[i,j] = 0
      }
    }
  }
  return(kRNG)
}

####################################################### Index to compare two adjacency matrices

GRI <- function(G1,G2){
  n = dim(G1)[1]
  return(1-(norm(G1-G2,type="F")^2/(n*(n-1))))
}

Reciprocal <-
  function(x){
    
    n = length(x)
    temp = c()
    for(i in 1: n){
      if(x[i] == 0) temp[i] = 0
      else temp[i] = 1/x[i]
    }
    return(temp)
  }

laplacian <-
  function(A, normalised = FALSE){
    
    n = dim(A)[1]
    temp = apply(abs(A), 2, sum)
    D = diag(temp, nrow = n)
    
    # temp1 = Reciprocal(sqrt(temp))
    temp1 = Reciprocal(temp)
    inv.D = diag(temp1, nrow = n)
    # half.D = diag(temp1, nrow = n)
    # if(normalised == TRUE) 	return(solve(D) %*% (D - A))
    if(normalised == TRUE) 	return(inv.D %*% (D - A))
    if(normalised == FALSE) return(D - A)
    
  }

spectral.clustering.rw <- function(A, normalised = TRUE, score = FALSE, K = 2, adj = FALSE){
  
  ### preparation 
  n = dim(A)[1]
  iso.A = isolate(A)
  iso.seq = iso.A$isolate
  noniso.seq = iso.A$nonisolate
  A.noniso = A[noniso.seq, noniso.seq]
  labels = rep(0, n)
  
  ### svd
  n = dim(A.noniso)[1]
  eig.A = eigen(A.noniso)
  if(score == F){
    U = matrix(0, nrow = n, ncol = K)	
  }
  if(score == T){
    U = matrix(0, nrow = n, ncol = (K-1))
  }
  
  ### get U matrix
  if(adj == F & score == F){
    L = laplacian(A = A.noniso, normalised = normalised)
    eig = eigen(L, symmetric = TRUE)
    
    for(i in 1:K){
      U[, i] = eig$vectors[, (n + 1 - i)]
    }
  }
  
  if(adj == T & score == F){
    ordered.vec = eig.A$vectors[, order(abs(eig.A$values), decreasing = T)]
    for(i in 1:K){
      U[, i] = ordered.vec[, (n + 1 - i)]
    }
  }
  
  if(score == T){
    ordered.vec = eig.A$vectors[, order(abs(eig.A$values), decreasing = T)]
    benchmark = ordered.vec[,1] + 1e-5
    for(i in 2:K){
      U[, (i-1)] = eig.A$vectors[, order(abs(eig.A$values), decreasing = F)][, i]/benchmark
    }
  }
  
  ### k-means
  #U = scale(U, center = F)
  temp = unique(U, margin = 2)
  if(dim(temp)[1] < K){stop('FAIL!')}
  k.means = kmeans(U, centers = K)
  labels[noniso.seq] = k.means$cluster
  
  return(labels)
}
