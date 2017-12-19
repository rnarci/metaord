import_data <- function(size_fraction, samples = NULL){
oldwd <- getwd()  
  
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

############################################################ Import data (note : pour l'instant, je ne prends pas les 
############################################################ matrices de distance asymetriques et celle de 
############################################################ sorensen-braycurtis PA)

if(size_fraction=="0.8-5" | size_fraction=="5-20" | size_fraction=="180-2000"){
  jaccard_abundance2 = read.table(file="mat_abundance_jaccard.csv",sep="")
  if(dim(jaccard_abundance2)[2]==1){
    jaccard_abundance2 = read.table(file="mat_abundance_jaccard.csv",sep=";")
    sep = ";"
  }else{
    sep = ""
  }
  delete=c(which(duplicated(jaccard_abundance2)==TRUE)) # remove duplicates
  jaccard_abundance2 = jaccard_abundance2[-delete,-delete]
  metagenomic_sample <<- as.character(jaccard_abundance2[-1,1])
  
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
  
  ############################################################ Distance matrices
  
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
  
  jaccard_abundance <<- matrix(NA,nrow=n,ncol=n)
  ab_jaccard_abundance <<- matrix(NA,nrow=n,ncol=n)
  braycurtis_abundance <<- matrix(NA,nrow=n,ncol=n)
  ab_ochiai_abundance <<- matrix(NA,nrow=n,ncol=n)
  ab_sorensen_abundance <<- matrix(NA,nrow=n,ncol=n)
  simka_jaccard_abundance <<- matrix(NA,nrow=n,ncol=n)
  
  chord_prevalence <<- matrix(NA,nrow=n,ncol=n)
  jaccard_prevalence <<- matrix(NA,nrow=n,ncol=n)
  kulczynski_prevalence <<- matrix(NA,nrow=n,ncol=n)
  ochiai_prevalence <<- matrix(NA,nrow=n,ncol=n)
  whittaker_prevalence <<- matrix(NA,nrow=n,ncol=n)
  simka_jaccard_prevalence <<- matrix(NA,nrow=n,ncol=n)
  # sorensen_braycurtis_prevalence = matrix(NA,nrow=n,ncol=n)
  
  for(j in 1:n){
    jaccard_abundance[,j] <<- as.numeric(jaccard_abundance2[,j])
    ab_jaccard_abundance[,j] <<- as.numeric(ab_jaccard_abundance2[,j])
    braycurtis_abundance[,j] <<- as.numeric(braycurtis_abundance2[,j])
    ab_ochiai_abundance[,j] <<- as.numeric(ab_ochiai_abundance2[,j])
    ab_sorensen_abundance[,j] <<- as.numeric(ab_sorensen_abundance2[,j])
    simka_jaccard_abundance[,j] <<- as.numeric(simka_jaccard_abundance2[,j])
    
    chord_prevalence[,j] <<- as.numeric(chord_prevalence2[,j])
    jaccard_prevalence[,j] <<- as.numeric(jaccard_prevalence2[,j])
    kulczynski_prevalence[,j] <<- as.numeric(kulczynski_prevalence2[,j])
    ochiai_prevalence[,j] <<- as.numeric(ochiai_prevalence2[,j])
    whittaker_prevalence[,j] <<- as.numeric(whittaker_prevalence2[,j])
    simka_jaccard_prevalence[,j] <<- as.numeric(simka_jaccard_prevalence2[,j])
    # sorensen_braycurtis_prevalence[,j] = as.numeric(sorensen_braycurtis_prevalence2[,j])
  }
  
  ############################################################ Name rows and columns correctly
  
  matrices.list <- c("jaccard_abundance", "ab_jaccard_abundance", "braycurtis_abundance", 
                     "ab_ochiai_abundance", "ab_sorensen_abundance", "simka_jaccard_abundance", 
                     "chord_prevalence", "jaccard_prevalence", "kulczynski_prevalence", 
                     "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence")
  for (mat in matrices.list) {
    dist.mat <- get(mat)
    rownames(dist.mat) <- colnames(dist.mat) <- metagenomic_sample
    if (!is.null(samples)) {
      complete_sample <- intersect(metagenomic_sample, samples)
      dist.mat <- dist.mat[complete_sample, complete_sample]
    }
    assign(x = mat, value = dist.mat, envir = .GlobalEnv)
  }
  
  if (!is.null(samples)) {
    metagenomic_sample <<- complete_sample
    if (length(metagenomic_sample) == 0) {
      warning("No sample selected, check the value of `samples` argument")
    }
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
  
  metagenomic_sample <<- row.names(jaccard_abundance)
  
  ############################################################ Distance matrices
  
  jaccard_abundance <<- unname(as.matrix(jaccard_abundance))
  ochiai_abundance <<- unname(as.matrix(ochiai_abundance))
  sorensen_abundance <<- unname(as.matrix(sorensen_abundance))
  simka_jaccard_abundance <<- unname(as.matrix(simka_jaccard_abundance))
  
  
  chord_hellinger_prevalence <<- unname(as.matrix(chord_hellinger_prevalence))
  jaccard_canberra_prevalence <<- unname(as.matrix(jaccard_canberra_prevalence))
  kulczynski_prevalence <<- unname(as.matrix(kulczynski_prevalence))
  ochiai_prevalence <<- unname(as.matrix(ochiai_prevalence))
  whittaker_prevalence <<- unname(as.matrix(whittaker_prevalence))
  simka_jaccard_prevalence <<- unname(as.matrix(simka_jaccard_prevalence))
  sorensen_braycurtis_prevalence <<- unname(as.matrix(sorensen_braycurtis_prevalence))
  
  ############################################################ Name rows and columns correctly
  
  matrices.list <- c("jaccard_abundance", "ochiai_abundance", "sorensen_abundance", 
                      "simka_jaccard_abundance", "chord_hellinger_prevalence", 
                     "jaccard_canberra_prevalence", "kulczynski_prevalence", 
                     "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence",
                     "sorensen_braycurtis_prevalence")
  
  for (mat in matrices.list) {
    dist.mat <- get(mat)
    rownames(dist.mat) <- colnames(dist.mat) <- metagenomic_sample
    if (!is.null(samples)) {
      complete_sample <- intersect(metagenomic_sample, samples)
      dist.mat <- dist.mat[complete_sample, complete_sample]
    }
    assign(x = mat, value = dist.mat, envir = .GlobalEnv)
  }
  
  if (!is.null(samples)) {
    metagenomic_sample <<- complete_sample
    if (length(metagenomic_sample) == 0) {
      warning("No sample selected, check the value of `samples` argument")
    }
  }
}

############################################################ Restore working directory

setwd(oldwd)

}


# ############################## Collect all statements from distance matrix (symetric)
# 
# collect_all_statements_distance <- function(L) {
#   n = dim(L)[1]
#   S = matrix(NA,nrow=n^3,ncol=4) # statements
#   statement = 1
#   
#   for(i in 1:(n-2)){
#       for (j in (i+1):(n-1)){
#         for(k in (j+1):n){
#           if(L[i,j]<L[j,k] & L[i,k]<L[j,k]){
#             S[statement,1] = i
#             S[statement,2:4] = sort(c(i,j,k))
#             statement = statement + 1           
#           }
#           if(L[j,i]<L[i,k] & L[j,k]<L[i,k]){
#             S[statement,1] = j
#             S[statement,2:4] = sort(c(i,j,k))
#             statement = statement + 1           
#           }
#           if(L[k,i]<L[i,j] & L[k,j]<L[i,j]){
#             S[statement,1] = k
#             S[statement,2:4] = sort(c(i,j,k))
#             statement = statement + 1           
#           }
#         }
#       }
#     }
#   return(na.omit(S))
#   }
# 
# 
# ############################## Create two-moons data set
# 
# .numObjects<-function(numObjects,numClusters){
#   if(length(numObjects)<numClusters){
#     numObjects<-rep(numObjects,numClusters)
#   }
#   if(length(numObjects)>numClusters){
#     numObjects<-numObjects[1:numClusters]
#   }
#   numObjects
# }
# 
# .toCsv<-function(file,data,klasy,outputRowNames,outputColNames,csv2){
#   if(paste(file,"",sep="")!=""){
#     if(csv2){
#       write.table(cbind(1:dim(data)[1],klasy,data),file=file,sep=";",dec=",",row.names=outputRowNames,col.names=outputColNames)
#     }
#     else{
#       write.table(cbind(1:dim(data)[1],klasy,data),file=file,row.names=outputRowNames,col.names=outputColNames)
#     }
#   }
# }
# 
# shapes.two.moon<-function(numObjects,shape1a,shape2b,shape1rFrom, shape1rTo,shape2rFrom, shape2rTo, outputCsv, outputCsv2, outputColNames, outputRowNames){
#   lo<-.numObjects(numObjects,2)
#   x <- matrix(0, nrow=sum(lo), ncol=2)
#   
#   for(i in 1:sum(lo)){
#     alpha<-runif(1,0,2*pi)
#     if(i>lo[1]){
#       r=runif(1,shape2rFrom,shape2rTo)
#     }
#     else{
#       r=runif(1,shape1rFrom,shape1rTo)
#     }
#     x[i,1]<-r*cos(alpha)
#     x[i,2]<-r*sin(alpha)
#     if(i<=lo[1]){
#       x[i,1]=shape1a+abs(x[i,1])
#     }
#     else{
#       x[i,1]=-abs(x[i,1])
#       x[i,2]=x[i,2]-shape2b
#     }
#   }
#   data<-x
#   klasy<-c(rep(1,lo[1]),rep(2,lo[2]))
#   .toCsv(outputCsv,data,klasy,outputColNames, outputRowNames,FALSE)
#   .toCsv(outputCsv2,data,klasy, outputColNames, outputRowNames,TRUE)
#   list(data=data,clusters=klasy)
# }
# 
# ############################## Algorithm 5 Clustering : unweighted version
# 
# k_RNG_clustering_unweighted <- function(S,k,l,n_data){
#   n_statements = dim(S)[1]
#   n_data = n_data
#   N = matrix(0,nrow=n_data,ncol=n_data)
#   D = matrix(0,nrow=n_data,ncol=n_data)
#   V = matrix(NA,nrow=n_data,ncol=n_data)
#   W = matrix(NA,nrow=n_data,ncol=n_data)
#   
#   for(i in 1:n_data){
#     for(j in 1:n_data){
#       for(y in 1:n_statements){
#         if(i!=j & (S[y,2]==i | S[y,3]==i | S[y,4]==i) & (S[y,2]==j | S[y,3]==j | S[y,4]==j)){
#           D[i,j] = D[i,j]+1
#           if(S[y,1]!=i & S[y,1]!=j){
#             N[i,j]=N[i,j]+1
#           }
#         }
#       }
#       if(D[i,j]==0){
#         V[i,j] = 1000000 
#       }else{
#         V[i,j] = N[i,j]/D[i,j]
#       }
#       if(V[i,j]<k/(n_data-2)){
#         W[i,j] = 1 
#       }else{
#         W[i,j] = 0
#       }
#     }
#   }
#   res = spectral.clustering(W, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   return(list(res=res,kRNG=W))
# }
# 
# 
# ############################## Algorithm 5 Clustering : weighted version
# 
# k_RNG_clustering_weighted <- function(S,k,l,n_data,sigma){
#   n_statements = dim(S)[1]
#   n_data = n_data
#   sigma = sigma
#   N = matrix(0,nrow=n_data,ncol=n_data)
#   D = matrix(0,nrow=n_data,ncol=n_data)
#   V = matrix(NA,nrow=n_data,ncol=n_data)
#   W = matrix(NA,nrow=n_data,ncol=n_data)
#   
#   for(i in 1:n_data){
#     for(j in 1:n_data){
#       for(y in 1:n_statements){
#         if(i!=j & (S[y,2]==i | S[y,3]==i | S[y,4]==i) & (S[y,2]==j | S[y,3]==j | S[y,4]==j)){
#           D[i,j] = D[i,j]+1
#           if(S[y,1]!=i & S[y,1]!=j){
#             N[i,j]=N[i,j]+1
#           }
#         }
#       }
#       if(D[i,j]==0){
#         V[i,j] = 1000000 
#       }else{
#         V[i,j] = N[i,j]/D[i,j]
#       }
#       if(V[i,j]<k/(n_data-2)){
#         W[i,j] = exp(-(V[i,j]^2)/(sigma^2)) 
#       }else{
#         W[i,j] = 0
#       }
#     }
#   }
#   res = spectral.clustering(W, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   return(list(res=res,kRNG=W))
# }
# 
# 
# ############################## Clustering with true k-RNG
# 
# # true_kRNG_clustering <- function(L,k,l) {
# #   n = dim(L)[1]
# #   count = 0
# #   kRNG = matrix(NA,nrow=n,ncol=n)
# # 
# #   for(i in 1:n){
# #     for (j in 1:n){
# #       for(m in 1:n){
# #         if(i!=j & i!=m & j!=m & L[i,j]>max(L[i,m],L[j,m])){
# #           count = count+1
# #         }
# #       }
# #       if(count<k & i!=j){
# #         kRNG[i,j] = 1
# #         count = 0
# #       }
# #       if(count>=k & i!=j){
# #         kRNG[i,j] = 0
# #         count = 0
# #       }
# #       if(i==j){
# #         kRNG[i,j] = 0
# #       }
# #     }
# #   }
# #   res = spectral.clustering(kRNG, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
# #   return(list(res=res,kRNG=kRNG))
# # }
# 
# true_kRNG_clustering_unweighted <- function(L,k,l) {
#   n = dim(L)[1]
#   count = 0
#   kRNG = matrix(NA,nrow=n,ncol=n)
#   
#   for(i in 1:n){
#     for (j in 1:n){
#       if(i!=j & length(which((L[i,j]>pmax(L[i,-c(i,j)],L[j,-c(i,j)]))==TRUE))<k){
#         kRNG[i,j] = 1
#       }
#       else{
#         kRNG[i,j] = 0
#       }
#     }
#   }
#   res = spectral.clustering(kRNG, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   return(list(res=res,kRNG=kRNG))
# }
# 
# true_kRNG_clustering_weighted <- function(L,k,l,sigma) {
#   n = dim(L)[1]
#   count = 0
#   kRNG = matrix(NA,nrow=n,ncol=n)
#   
#   for(i in 1:n){
#     for (j in 1:n){
#       N = length(which((L[i,j]>pmax(L[i,-c(i,j)],L[j,-c(i,j)]))==TRUE))
#       if(i!=j & N<k){
#         kRNG[i,j] = exp(-N^2/((n-2)^2*sigma^2))
#       }
#       else{
#         kRNG[i,j] = 0
#       }
#     }
#   }
#   res = spectral.clustering(kRNG, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
#   return(list(res=res,kRNG=kRNG))
# }
# 
# make.kRNG <- function(L,k,sigma){
#   n = dim(L)[1]
#   count = 0
#   kRNG = matrix(NA,nrow=n,ncol=n)
#   
#   for(i in 1:n){
#     for (j in 1:n){
#       N = length(which((L[i,j]>pmax(L[i,-c(i,j)],L[j,-c(i,j)]))==TRUE))
#       if(i!=j & N<k){
#         kRNG[i,j] = 1 #exp(-N^2/((n-2)^2*sigma^2))
#       }
#       else{
#         kRNG[i,j] = 0
#       }
#     }
#   }
#   return(kRNG)
# }
# 
# ####################################################### Index to compare two adjacency matrices
# 
# GRI <- function(G1,G2){
#   n = dim(G1)[1]
#   return(1-(norm(G1-G2,type="F")^2/(n*(n-1))))
# }
# 
# ####################################################### Spectral clustering
# 
# Reciprocal <-
#   function(x){
#     
#     n = length(x)
#     temp = c()
#     for(i in 1: n){
#       if(x[i] == 0) temp[i] = 0
#       else temp[i] = 1/x[i]
#     }
#     return(temp)
#   }
# 
# laplacian <-
#   function(A, normalised = FALSE){
#     
#     n = dim(A)[1]
#     temp = apply(abs(A), 2, sum)
#     D = diag(temp, nrow = n)
#     
#     # temp1 = Reciprocal(sqrt(temp))
#     temp1 = Reciprocal(temp)
#     inv.D = diag(temp1, nrow = n)
#     # half.D = diag(temp1, nrow = n)
#     # if(normalised == TRUE) 	return(solve(D) %*% (D - A))
#     if(normalised == TRUE) 	return(inv.D %*% (D - A))
#     if(normalised == FALSE) return(D - A)
#     
#   }
# 
# spectral.clustering.rw <- function(A, normalised = TRUE, score = FALSE, K = 2, adj = FALSE){
#   
#   ### preparation 
#   n = dim(A)[1]
#   iso.A = isolate(A)
#   iso.seq = iso.A$isolate
#   noniso.seq = iso.A$nonisolate
#   A.noniso = A[noniso.seq, noniso.seq]
#   labels = rep(0, n)
#   
#   ### svd
#   n = dim(A.noniso)[1]
#   eig.A = eigen(A.noniso)
#   if(score == F){
#     U = matrix(0, nrow = n, ncol = K)	
#   }
#   if(score == T){
#     U = matrix(0, nrow = n, ncol = (K-1))
#   }
#   
#   ### get U matrix
#   if(adj == F & score == F){
#     L = laplacian(A = A.noniso, normalised = normalised)
#     eig = eigen(L, symmetric = TRUE)
#     
#     for(i in 1:K){
#       U[, i] = eig$vectors[, (n + 1 - i)]
#     }
#   }
#   
#   if(adj == T & score == F){
#     ordered.vec = eig.A$vectors[, order(abs(eig.A$values), decreasing = T)]
#     for(i in 1:K){
#       U[, i] = ordered.vec[, (n + 1 - i)]
#     }
#   }
#   
#   if(score == T){
#     ordered.vec = eig.A$vectors[, order(abs(eig.A$values), decreasing = T)]
#     benchmark = ordered.vec[,1] + 1e-5
#     for(i in 2:K){
#       U[, (i-1)] = eig.A$vectors[, order(abs(eig.A$values), decreasing = F)][, i]/benchmark
#     }
#   }
#   
#   ### k-means
#   #U = scale(U, center = F)
#   temp = unique(U, margin = 2)
#   if(dim(temp)[1] < K){stop('FAIL!')}
#   k.means = kmeans(U, centers = K)
#   labels[noniso.seq] = k.means$cluster
#   
#   return(labels)
# }