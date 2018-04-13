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

############################################################ Import data 

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
  sorensen_braycurtis_prevalence2 = read.table(file="mat_presenceAbsence_sorensen-braycurtis.csv",sep="\t") 
  
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
  sorensen_braycurtis_prevalence2 = unname(as.matrix(sorensen_braycurtis_prevalence2)[c(-1,-delete),c(-1,-delete)])
  
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
  sorensen_braycurtis_prevalence <<- matrix(NA,nrow=n,ncol=n)
  
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
    sorensen_braycurtis_prevalence[,j] <<- as.numeric(sorensen_braycurtis_prevalence2[,j])
  }
  
  ############################################################ Name rows and columns correctly
  
  matrices.list <<- c("jaccard_abundance", "ab_jaccard_abundance", "braycurtis_abundance", 
                     "ab_ochiai_abundance", "ab_sorensen_abundance", "simka_jaccard_abundance", 
                     "chord_prevalence", "jaccard_prevalence", "kulczynski_prevalence", 
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
  
}else{
  sep = ""
  jaccard_abundance2 = read.table(file="mat_abundance_jaccard.csv",sep=sep)
  ochiai_abundance2 = read.table(file="mat_abundance_ochiai.csv",sep=sep)
  sorensen_abundance2 = read.table(file="mat_abundance_sorensen.csv",sep=sep)
  simka_jaccard_abundance2 = read.table(file="mat_abundance_simka-jaccard.csv",sep=sep)
  
  chord_hellinger_prevalence2 = read.table(file="mat_presenceAbsence_chord-hellinger.csv",sep=sep)
  jaccard_canberra_prevalence2 = read.table(file="mat_presenceAbsence_jaccard-canberra.csv",sep=sep)
  kulczynski_prevalence2 = read.table(file="mat_presenceAbsence_kulczynski.csv",sep=sep)
  ochiai_prevalence2 = read.table(file="mat_presenceAbsence_ochiai.csv",sep=sep)
  whittaker_prevalence2 = read.table(file="mat_presenceAbsence_whittaker.csv",sep=sep)
  simka_jaccard_prevalence2 = read.table(file="mat_presenceAbsence_simka-jaccard.csv",sep=sep)
  sorensen_braycurtis_prevalence2 = read.table(file="mat_presenceAbsence_sorensen-braycurtis.csv",sep=sep)
  
  metagenomic_sample <<- row.names(jaccard_abundance2)
  
  ############################################################ Distance matrices
  
  jaccard_abundance2 = unname(as.matrix(jaccard_abundance2))
  ochiai_abundance2 = unname(as.matrix(ochiai_abundance2))
  sorensen_abundance2 = unname(as.matrix(sorensen_abundance2))
  simka_jaccard_abundance2 = unname(as.matrix(simka_jaccard_abundance2))
  
  
  chord_hellinger_prevalence2 = unname(as.matrix(chord_hellinger_prevalence2))
  jaccard_canberra_prevalence2 = unname(as.matrix(jaccard_canberra_prevalence2))
  kulczynski_prevalence2 = unname(as.matrix(kulczynski_prevalence2))
  ochiai_prevalence2 = unname(as.matrix(ochiai_prevalence2))
  whittaker_prevalence2 = unname(as.matrix(whittaker_prevalence2))
  simka_jaccard_prevalence2 = unname(as.matrix(simka_jaccard_prevalence2))
  sorensen_braycurtis_prevalence2 = unname(as.matrix(sorensen_braycurtis_prevalence2))
  
  n = dim(jaccard_abundance2)[1]
  
  jaccard_abundance <<- matrix(NA,nrow=n,ncol=n)
  ochiai_abundance <<- matrix(NA,nrow=n,ncol=n)
  sorensen_abundance <<- matrix(NA,nrow=n,ncol=n)
  simka_jaccard_abundance <<- matrix(NA,nrow=n,ncol=n)
  
  
  chord_hellinger_prevalence <<- matrix(NA,nrow=n,ncol=n)
  jaccard_canberra_prevalence <<- matrix(NA,nrow=n,ncol=n)
  kulczynski_prevalence <<- matrix(NA,nrow=n,ncol=n)
  ochiai_prevalence <<- matrix(NA,nrow=n,ncol=n)
  whittaker_prevalence <<- matrix(NA,nrow=n,ncol=n)
  simka_jaccard_prevalence <<- matrix(NA,nrow=n,ncol=n)
  sorensen_braycurtis_prevalence <<- matrix(NA,nrow=n,ncol=n)
  
  for(j in 1:n){
    jaccard_abundance[,j] <<- as.numeric(jaccard_abundance2[,j])
    ochiai_abundance[,j] <<- as.numeric(ochiai_abundance2[,j])
    sorensen_abundance[,j] <<- as.numeric(sorensen_abundance2[,j])
    simka_jaccard_abundance[,j] <<- as.numeric(simka_jaccard_abundance2[,j])
    
    chord_hellinger_prevalence[,j] <<- as.numeric(chord_hellinger_prevalence2[,j])
    jaccard_canberra_prevalence[,j] <<- as.numeric(jaccard_canberra_prevalence2[,j])
    kulczynski_prevalence[,j] <<- as.numeric(kulczynski_prevalence2[,j])
    ochiai_prevalence[,j] <<- as.numeric(ochiai_prevalence2[,j])
    whittaker_prevalence[,j] <<- as.numeric(whittaker_prevalence2[,j])
    simka_jaccard_prevalence[,j] <<- as.numeric(simka_jaccard_prevalence2[,j])
    sorensen_braycurtis_prevalence[,j] <<- as.numeric(sorensen_braycurtis_prevalence2[,j])
  }
  
  ############################################################ Name rows and columns correctly
  
  matrices.list <<- c("jaccard_abundance", "ochiai_abundance", "sorensen_abundance", 
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

get_data_RI <- function(){
  n = dim(jaccard_abundance)[1]
  nb_mat = length(matrices.list)
  
  ####################################### To get all of the clusterings by an ordinal method (k-NNG + spectral clustering)
  ####################################### for each distances matrix and each size fraction
  
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
  kNN13 = matrix(NA, nrow = n, ncol = n) 
  
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
  res13 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
  
  
  kNN = list(kNN1,kNN2,kNN3,kNN4,kNN5,kNN6,kNN7,kNN8,kNN9,kNN10,kNN11,kNN12,kNN13)
  res = list(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12,res13)
  
  a = 0
  
  for(k in 1:(round(n/2))){
    for(num_matrix in 1:nb_mat){
      kNN[[num_matrix]] <- make.kNNG(get(matrices.list[num_matrix]), k = k, symm = TRUE, weight = FALSE)
      for(l in 2:15){
        res[[num_matrix]][14*(k-1)+l-1,1] = k 
        res[[num_matrix]][14*(k-1)+l-1,2] = l
        res[[num_matrix]][14*(k-1)+l-1,3:(n+2)] = spectral.clustering.new(kNN[[num_matrix]], normalised = TRUE, score = FALSE, K = l, adj = FALSE)
      }
    }
    a = a + 14
  }
  
  for(num_matrix in 1:nb_mat){
    res[[num_matrix]] = as.data.frame(res[[num_matrix]])
  }
  
  if(size_fraction=="0.8-5"){
    setwd(dir="~/metaord/0.8-5/mat_spectral") 
  }
  
  if(size_fraction=="5-20"){
    setwd(dir="~/metaord/5-20/mat_spectral") 
  }
  
  if(size_fraction=="180-2000"){
    setwd(dir="~/metaord/180-2000/mat_spectral") 
  }
  
  if(size_fraction=="0.22-3"){
    setwd(dir="~/metaord/0.22-3/mat_spectral") 
  }
  
  if(size_fraction=="0-0.2"){
    setwd(dir="~/metaord/0-0.2/mat_spectral") 
  }
  
  if(size_fraction=="20-180"){
    setwd(dir="~/metaord/20-180/mat_spectral") 
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
  if(nb_mat == 13){
    write.table(x = res[[12]], file = "Mat12_spectral") 
    write.table(x = res[[13]], file = "Mat13_spectral") 
  }
  
  ####################################### To get all of the clusterings by a classical method (MDS + K-means)
  ####################################### for each distances matrix and each size fraction
  
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
  res13 = matrix(NA, nrow = round(n/2)*length(2:15), ncol = n + 2)
  
  res = list(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12,res13)
  
  a = 0
  
  for(k in 1:(round(n/2))){
    if(nb_mat == 13){
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
      fit13 = cmdscale(sorensen_braycurtis_prevalence, k=k)
    }else{
      fit1 = cmdscale(jaccard_abundance, k=k)
      fit2 = cmdscale(ochiai_abundance, k=k)
      fit3 = cmdscale(sorensen_abundance, k=k)
      fit4 = cmdscale(simka_jaccard_abundance, k=k)
      fit5 = cmdscale(chord_hellinger_prevalence, k=k)
      fit6 = cmdscale(jaccard_canberra_prevalence, k=k)
      fit7 = cmdscale(kulczynski_prevalence, k=k)
      fit8 = cmdscale(ochiai_prevalence, k=k)
      fit9 = cmdscale(whittaker_prevalence, k=k)
      fit10 = cmdscale(simka_jaccard_prevalence, k=k)
      fit11 = cmdscale(sorensen_braycurtis_prevalence, k=k)
    }
    for(l in 2:15){
      for(num_matrix in 1:nb_mat){
        res[[num_matrix]][14*(k-1)+l-1,1] = k
        res[[num_matrix]][14*(k-1)+l-1,2] = l
        res[[num_matrix]][14*(k-1)+l-1,3:(n+2)] = kmeans(x = get(paste("fit",num_matrix,sep="")), centers = l, nstart = 1000, iter.max = 20)$cluster
      }
    }
    a = a + 14
  }
  
  for(num_matrix in 1:nb_mat){
    res[[num_matrix]] = as.data.frame(res[[num_matrix]])
  }

  if(size_fraction=="0.8-5"){
    setwd(dir="~/metaord/0.8-5/mat_kmeans")
  }

  if(size_fraction=="5-20"){
    setwd(dir="~/metaord/5-20/mat_kmeans")
  }

  if(size_fraction=="180-2000"){
    setwd(dir="~/metaord/180-2000/mat_kmeans")
  }

  if(size_fraction=="0.22-3"){
    setwd(dir="~/metaord/0.22-3/mat_kmeans")
  }

  if(size_fraction=="0-0.2"){
    setwd(dir="~/metaord/0-0.2/mat_kmeans")
  }

  if(size_fraction=="20-180"){
    setwd(dir="~/metaord/20-180/mat_kmeans")
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
  if(nb_mat == 13){
    write.table(x = res[[12]], file = "Mat12_kmeans")
    write.table(x = res[[13]], file = "Mat13_kmeans")
  }

  ####################################### To get all of the Rand indices between clusterings for a given number of clusters and for 
  ####################################### each pair of values of k (number of nearest neighbours)
  
  if(size_fraction=="0.8-5"){
    setwd(dir="~/metaord/0.8-5/mat_spectral")
  }

  if(size_fraction=="5-20"){
    setwd(dir="~/metaord/5-20/mat_spectral")
  }

  if(size_fraction=="180-2000"){
    setwd(dir="~/metaord/180-2000/mat_spectral")
  }

  if(size_fraction=="0.22-3"){
    setwd(dir="~/metaord/0.22-3/mat_spectral")
  }

  if(size_fraction=="0-0.2"){
    setwd(dir="~/metaord/0-0.2/mat_spectral")
  }

  if(size_fraction=="20-180"){
    setwd(dir="~/metaord/20-180/mat_spectral")
  }

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
  if(nb_mat == 13){
    mat12 = read.table(file = "Mat12_spectral")
    mat13 = read.table(file = "Mat13_spectral")
    mat = list(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8,mat9,mat10,mat11,mat12,mat13)
  }else{
    mat = list(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8,mat9,mat10,mat11)
  }
  
  for(num_matrix in 1:nb_mat){
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
  if(nb_mat == 13){
    RI12 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
    RI13 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
    RI = list(RI1,RI2,RI3,RI4,RI5,RI6,RI7,RI8,RI9,RI10,RI11,RI12,RI13)
  }else{
    RI = list(RI1,RI2,RI3,RI4,RI5,RI6,RI7,RI8,RI9,RI10,RI11)
  }
  
  for(num_matrix in 1:nb_mat){
    a = 1
    for(k1 in 1:(round(n/2)-1)){
      for(k2 in (k1+1):round(n/2)){
        for(l in 2:15){
          RI[[num_matrix]][a,1] = l
          RI[[num_matrix]][a,2] = k1
          RI[[num_matrix]][a,3] = k2
          RI[[num_matrix]][a,4] = rand.index(mat[[num_matrix]][14*(k1-1)+l-1,3:(n+2)], mat[[num_matrix]][14*(k2-1)+l-1,3:(n+2)])
          a = a + 1
        }
      }
    }
  }
  
  
  for(num_matrix in 1:nb_mat){
    RI[[num_matrix]] = as.data.frame(RI[[num_matrix]])
  }
  
  res <- lapply(1:nb_mat, function(i) {
    res <- RI[[i]]
    colnames(res) <- c("Nb_Cluster", "k1", "k2", "RI")
    res %>% mutate(Matrice.Nb = i)
  }) %>%
    bind_rows() %>%
    mutate(Distance = matrices.list[Matrice.Nb])
  
  if(size_fraction=="0.8-5"){
    setwd(dir="~/metaord/0.8-5/RI_results")
  }

  if(size_fraction=="5-20"){
    setwd(dir="~/metaord/5-20/RI_results")
  }

  if(size_fraction=="180-2000"){
    setwd(dir="~/metaord/180-2000/RI_results")
  }

  if(size_fraction=="0.22-3"){
    setwd(dir="~/metaord/0.22-3/RI_results")
  }

  if(size_fraction=="0-0.2"){
    setwd(dir="~/metaord/0-0.2/RI_results")
  }

  if(size_fraction=="20-180"){
    setwd(dir="~/metaord/20-180/RI_results")
  }

  fwrite(res, file = "RI_spectral_all")

  ####################################### To get all of the Rand indices between clusterings for a given number of clusters and for 
  ####################################### each pair of values of d (dimension of the new space of the projected points)
  
  if(size_fraction=="0.8-5"){
    setwd(dir="~/metaord/0.8-5/mat_kmeans")
  }
  
  if(size_fraction=="5-20"){
    setwd(dir="~/metaord/5-20/mat_kmeans")
  }
  
  if(size_fraction=="180-2000"){
    setwd(dir="~/metaord/180-2000/mat_kmeans")
  }
  
  if(size_fraction=="0.22-3"){
    setwd(dir="~/metaord/0.22-3/mat_kmeans")
  }
  
  if(size_fraction=="0-0.2"){
    setwd(dir="~/metaord/0-0.2/mat_kmeans")
  }
  
  if(size_fraction=="20-180"){
    setwd(dir="~/metaord/20-180/mat_kmeans")
  }

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
  if(nb_mat == 13){
    mat12 = read.table(file = "Mat12_kmeans")
    mat13 = read.table(file = "Mat13_kmeans")
    mat = list(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8,mat9,mat10,mat11,mat12,mat13)
  }else{
    mat = list(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8,mat9,mat10,mat11)
  }
  
  for(num_matrix in 1:nb_mat){
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
  if(nb_mat == 13){
    RI12 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
    RI13 = matrix(NA, nrow = 14*round(n/2)*(round(n/2)-1)/2, ncol = 4)
    RI = list(RI1,RI2,RI3,RI4,RI5,RI6,RI7,RI8,RI9,RI10,RI11,RI12,RI13)
  }else{
    RI = list(RI1,RI2,RI3,RI4,RI5,RI6,RI7,RI8,RI9,RI10,RI11)
  }
  
  for(num_matrix in 1:nb_mat){
    a = 1
    for(k1 in 1:(round(n/2)-1)){
      for(k2 in (k1+1):round(n/2)){
        for(l in 2:15){
          RI[[num_matrix]][a,1] = l
          RI[[num_matrix]][a,2] = k1
          RI[[num_matrix]][a,3] = k2
          RI[[num_matrix]][a,4] = rand.index(mat[[num_matrix]][14*(k1-1)+l-1,3:(n+2)], mat[[num_matrix]][14*(k2-1)+l-1,3:(n+2)])
          a = a + 1
        }
      }
    }
  }
  
  
  for(num_matrix in 1:nb_mat){
    RI[[num_matrix]] = as.data.frame(RI[[num_matrix]])
  }
  
  res <- lapply(1:nb_mat, function(i) {
    res <- RI[[i]]
    colnames(res) <- c("Nb_Cluster", "d1", "d2", "RI")
    res %>% mutate(Matrice.Nb = i)
  }) %>%
    bind_rows() %>%
    mutate(Distance = matrices.list[Matrice.Nb])

  if(size_fraction=="0.8-5"){
    setwd(dir="~/metaord/0.8-5/RI_results")
  }

  if(size_fraction=="5-20"){
    setwd(dir="~/metaord/5-20/RI_results")
  }

  if(size_fraction=="180-2000"){
    setwd(dir="~/metaord/180-2000/RI_results")
  }

  if(size_fraction=="0.22-3"){
    setwd(dir="~/metaord/0.22-3/RI_results")
  }

  if(size_fraction=="0-0.2"){
    setwd(dir="~/metaord/0-0.2/RI_results")
  }

  if(size_fraction=="20-180"){
    setwd(dir="~/metaord/20-180/RI_results")
  }

  fwrite(res, file = "RI_kmeans_all")
  
}

compare_clusterings <- function(){
  
  ########################### Function comparing distinct clusterings by the calculation of the Rand index
  
  library(gtools)
  
  type <- function(i){
    if(even(i)){
      type = "Spectral"
    }else{
      type = "Kmeans"
    }
    return(type)
  }
  
  matrices.list1 = c("jaccard_abundance", "ab_jaccard_abundance", "braycurtis_abundance", 
                     "ab_ochiai_abundance", "ab_sorensen_abundance", "simka_jaccard_abundance", 
                     "chord_prevalence", "jaccard_prevalence", "kulczynski_prevalence", 
                     "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence",
                     "sorensen_braycurtis_prevalence")
  matrices.list2 = c("jaccard_abundance", "ochiai_abundance", "sorensen_abundance", 
                     "simka_jaccard_abundance", "chord_hellinger_prevalence", 
                     "jaccard_canberra_prevalence", "kulczynski_prevalence", 
                     "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence",
                     "sorensen_braycurtis_prevalence")
  matrices.list.spectral1 = c("Mat1_spectral", "Mat2_spectral", "Mat3_spectral", "Mat4_spectral", 
                              "Mat5_spectral", "Mat6_spectral", "Mat7_spectral", "Mat8_spectral", 
                              "Mat9_spectral", "Mat10_spectral", "Mat11_spectral", "Mat12_spectral", 
                              "Mat13_spectral")
  matrices.list.kmeans1 = c("Mat1_kmeans", "Mat2_kmeans", "Mat3_kmeans", "Mat4_kmeans", 
                            "Mat5_kmeans", "Mat6_kmeans", "Mat7_kmeans", "Mat8_kmeans", 
                            "Mat9_kmeans", "Mat10_kmeans", "Mat11_kmeans", "Mat12_kmeans", 
                            "Mat13_kmeans")
  
  matrices.list.spectral2  = c("Mat1_spectral", "Mat2_spectral", "Mat3_spectral", "Mat4_spectral", 
                               "Mat5_spectral", "Mat6_spectral", "Mat7_spectral", "Mat8_spectral", 
                               "Mat9_spectral", "Mat10_spectral", "Mat11_spectral")
  matrices.list.kmeans2    = c("Mat1_kmeans", "Mat2_kmeans", "Mat3_kmeans", "Mat4_kmeans", 
                               "Mat5_kmeans", "Mat6_kmeans", "Mat7_kmeans", "Mat8_kmeans", 
                               "Mat9_kmeans", "Mat10_kmeans", "Mat11_kmeans")
  
  size_fraction = "0.8-5"
  delete = which(is.na(data.0.8_5$Genocenose))
  data.0.8_5 <<- data.0.8_5[-delete,]
  import_data(size_fraction, samples = rownames(data.0.8_5))
  data.0.8_5 <<- data.0.8_5[metagenomic_sample,]
  
  size_fraction = "0-0.2"
  delete = which(is.na(data.0_0.2$Genocenose))
  data.0_0.2 <<- data.0_0.2[-delete,]
  import_data(size_fraction, samples = rownames(data.0_0.2))
  data.0_0.2 <<- data.0_0.2[metagenomic_sample,]
  
  size_fraction = "0.22-3"
  delete = which(is.na(data.0.22_3$Genocenose))
  data.0.22_3 <<- data.0.22_3[-delete,]
  import_data(size_fraction, samples = rownames(data.0.22_3))
  data.0.22_3 <<- data.0.22_3[metagenomic_sample,]
  
  size_fraction = "5-20"
  delete = which(is.na(data.5_20$Genocenose))
  data.5_20 <<- data.5_20[-delete,]
  import_data(size_fraction, samples = rownames(data.5_20))
  data.5_20 <<- data.5_20[metagenomic_sample,]
  
  size_fraction = "20-180"
  delete = which(is.na(data.20_180$Genocenose))
  data.20_180 <<- data.20_180[-delete,]
  import_data(size_fraction, samples = rownames(data.20_180))
  data.20_180 <<- data.20_180[metagenomic_sample,]
  
  size_fraction = "180-2000"
  # delete = which(is.na(data.180_2000$Genocenose))
  # data.180_2000 = data.180_2000[-delete,]
  import_data(size_fraction, samples = rownames(data.180_2000))
  data.180_2000 <<- data.180_2000[metagenomic_sample,]
  
  res_1 = data.frame()
  res_2 = data.frame()
  res_3 = data.frame()
  res_4 = data.frame()
  res_5 = data.frame()
  res_6 = data.frame()
  res_7 = data.frame()
  
  res1 = list(res_1,res_2,res_3,res_4,res_5,res_6,res_7)
  res2 = list(res_1,res_2,res_3,res_4,res_5,res_6)
  res3 = list(res_1,res_2,res_3,res_4)
  
  res = list(res1,res2,res3)
  
  a = 1
  
  ################################### Compare the clusterings obtained by Genoscope and those obtained 
  ################################### by the two methods (classical and ordinal) for distances matrix 
  ################################### which are present in all of the 6 sizes fraction
  
  for(select_distances_matrix in c("jaccard_abundance", "simka_jaccard_abundance", 
                                   "kulczynski_prevalence", "ochiai_prevalence",
                                   "whittaker_prevalence", "simka_jaccard_prevalence", 
                                   "sorensen_braycurtis_prevalence")){
    num_matrix1 = which(matrices.list1 == select_distances_matrix)
    num_matrix2 = which(matrices.list2 == select_distances_matrix)
    
    l1 = max(data.0.8_5$Genocenose)
    l2 = max(data.0_0.2$Genocenose)
    l3 = max(data.0.22_3$Genocenose)
    l4 = max(data.5_20$Genocenose)
    l5 = max(data.20_180$Genocenose)
    l6 = max(data.180_2000$Genocenose)
    
    n1 = dim(data.0.8_5)[1]
    n2 = dim(data.0_0.2)[1]
    n3 = dim(data.0.22_3)[1]
    n4 = dim(data.5_20)[1]
    n5 = dim(data.20_180)[1]
    n6 = dim(data.180_2000)[1]
    
    dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.8-5/mat_kmeans",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.8-5/mat_spectral",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0-0.2/mat_kmeans",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0-0.2/mat_spectral",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.22-3/mat_kmeans",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.22-3/mat_spectral",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/5-20/mat_kmeans",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/5-20/mat_spectral",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/20-180/mat_kmeans",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/20-180/mat_spectral",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/180-2000/mat_kmeans",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/180-2000/mat_spectral")
    
    fraction = c("0.8-5","0-0.2","0.22-3","5-20","20-180","180-2000")
    data = list(data.0.8_5$Genocenose,data.0_0.2$Genocenose,data.0.22_3$Genocenose,data.5_20$Genocenose,
                data.20_180$Genocenose,data.180_2000$Genocenose)
    l = list(l1,l2,l3,l4,l5,l6)
    n = list(n1,n2,n3,n4,n5,n6)
    
    RI1 = matrix(0, nrow = round(n[[1]]/2), ncol = 2)
    RI2 = matrix(0, nrow = round(n[[1]]/2), ncol = 2)
    RI3 = matrix(0, nrow = round(n[[2]]/2), ncol = 2)
    RI4 = matrix(0, nrow = round(n[[2]]/2), ncol = 2)
    RI5 = matrix(0, nrow = round(n[[3]]/2), ncol = 2)
    RI6 = matrix(0, nrow = round(n[[3]]/2), ncol = 2)
    RI7 = matrix(0, nrow = round(n[[4]]/2), ncol = 2)
    RI8 = matrix(0, nrow = round(n[[4]]/2), ncol = 2)
    RI9 = matrix(0, nrow = round(n[[5]]/2), ncol = 2)
    RI10 = matrix(0, nrow = round(n[[5]]/2), ncol = 2)
    RI11 = matrix(0, nrow = round(n[[6]]/2), ncol = 2)
    RI12 = matrix(0, nrow = round(n[[6]]/2), ncol = 2)
    
    RI = list(RI1,RI2,RI3,RI4,RI5,RI6,RI7,RI8,RI9,RI10,RI11,RI12)
    
    for(num_RI in 1:12){
      setwd(dir=dir[num_RI])
      if(even(num_RI) & num_RI == 2 | num_RI == 8 |num_RI == 12){
        mat_spectral = read.table(file = matrices.list.spectral1[num_matrix1])
        mat = as.matrix(mat_spectral) 
      }
      if(even(num_RI) & num_RI == 4 | num_RI == 6 |num_RI == 10){
        mat_spectral = read.table(file = matrices.list.spectral2[num_matrix2])
        mat = as.matrix(mat_spectral) 
      }
      if(odd(num_RI) & num_RI == 1 | num_RI == 7 |num_RI == 11){
        mat_kmeans = read.table(file = matrices.list.kmeans1[num_matrix1])
        mat = as.matrix(mat_kmeans) 
      }
      if(odd(num_RI) & num_RI == 3 | num_RI == 5 |num_RI == 9){
        mat_kmeans = read.table(file = matrices.list.kmeans2[num_matrix2])
        mat = as.matrix(mat_kmeans) 
      }
      for(k in 1:round(n[[ceiling(num_RI/2)]]/2)){
        RI[[num_RI]][k,1] = k
        RI[[num_RI]][k,2] = rand.index(mat[14*(k-1)+l[[ceiling(num_RI/2)]]-1,3:(n[[ceiling(num_RI/2)]]+2)], data[[ceiling(num_RI/2)]])
      } 
      RI[[num_RI]] = as.data.frame(RI[[num_RI]])
    }
    
    res1[[a]] <- lapply(1:12, function(i) {
      res1[[a]] <- RI[[i]]
      colnames(res1[[a]]) <- c("Param", "RI")
      res1[[a]] %>% mutate(Type = type(i), Fraction = fraction[ceiling(i/2)], Distance = select_distances_matrix)
    }) %>%
      bind_rows() 
    a = a + 1
  }
  
  res1 = bind_rows(res1[[1]], res1[[2]], res1[[3]], res1[[4]], res1[[5]], res1[[6]], res1[[7]])
  a = 1
  
  ################################### Compare the clusterings obtained by Genoscope and those obtained 
  ################################### by the two methods (classical and ordinal) for distances matrix 
  ################################### which are only present in 3 of the 6 sizes fraction
  
  for(select_distances_matrix in c("ab_jaccard_abundance", "braycurtis_abundance",
                                   "ab_ochiai_abundance", "ab_sorensen_abundance",
                                   "chord_prevalence", "jaccard_prevalence")){
    num_matrix1 = which(matrices.list1 == select_distances_matrix)
    
    l1 = max(data.0.8_5$Genocenose)
    l2 = max(data.5_20$Genocenose)
    l3 = max(data.180_2000$Genocenose)
    
    n1 = dim(data.0.8_5)[1]
    n2 = dim(data.5_20)[1]
    n3 = dim(data.180_2000)[1]
    
    dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.8-5/mat_kmeans",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.8-5/mat_spectral",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/5-20/mat_kmeans",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/5-20/mat_spectral",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/180-2000/mat_kmeans",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/180-2000/mat_spectral")
    fraction = c("0.8-5","5-20","180-2000")
    data = list(data.0.8_5$Genocenose,data.5_20$Genocenose,data.180_2000$Genocenose)
    l = list(l1,l2,l3)
    n = list(n1,n2,n3)
    
    RI1 = matrix(0, nrow = round(n[[1]]/2), ncol = 2)
    RI2 = matrix(0, nrow = round(n[[1]]/2), ncol = 2)
    RI3 = matrix(0, nrow = round(n[[2]]/2), ncol = 2)
    RI4 = matrix(0, nrow = round(n[[2]]/2), ncol = 2)
    RI5 = matrix(0, nrow = round(n[[3]]/2), ncol = 2)
    RI6 = matrix(0, nrow = round(n[[3]]/2), ncol = 2)
    
    RI = list(RI1,RI2,RI3,RI4,RI5,RI6)
    
    for(num_RI in 1:6){
      setwd(dir=dir[num_RI])
      if(even(num_RI)){
        mat_spectral = read.table(file = matrices.list.spectral1[num_matrix1])
        mat = as.matrix(mat_spectral) 
      }else{
        mat_kmeans = read.table(file = matrices.list.kmeans1[num_matrix1])
        mat = as.matrix(mat_kmeans) 
      }
      for(k in 1:round(n[[ceiling(num_RI/2)]]/2)){
        RI[[num_RI]][k,1] = k
        RI[[num_RI]][k,2] = rand.index(mat[14*(k-1)+l[[ceiling(num_RI/2)]]-1,3:(n[[ceiling(num_RI/2)]]+2)], data[[ceiling(num_RI/2)]])
      } 
      RI[[num_RI]] = as.data.frame(RI[[num_RI]])
    }
    
    res2[[a]] <- lapply(1:6, function(i) {
      res2[[a]] <- RI[[i]]
      colnames(res2[[a]]) <- c("Param", "RI")
      res2[[a]] %>% mutate(Type = type(i), Fraction = fraction[ceiling(i/2)], Distance = select_distances_matrix)
    }) %>%
      bind_rows() 
    a = a + 1
  }
  
  res2 = bind_rows(res2[[1]], res2[[2]], res2[[3]], res2[[4]], res2[[5]], res2[[6]])
  a = 1
  
  ################################### Compare the clusterings obtained by Genoscope and those obtained 
  ################################### by the two methods (classical and ordinal) for distances matrix 
  ################################### which are only present in the 3 remaining sizes fraction
  
  for(select_distances_matrix in c("ochiai_abundance", "sorensen_abundance", "chord_hellinger_prevalence",
                                   "jaccard_canberra_prevalence")){
    num_matrix2 = which(matrices.list2 == select_distances_matrix)

    l1 = max(data.0_0.2$Genocenose)
    l2 = max(data.0.22_3$Genocenose)
    l3 = max(data.20_180$Genocenose)
    
    n1 = dim(data.0_0.2)[1]
    n2 = dim(data.0.22_3)[1]
    n3 = dim(data.20_180)[1]
    
    dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0-0.2/mat_kmeans",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0-0.2/mat_spectral",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.22-3/mat_kmeans",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.22-3/mat_spectral",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/20-180/mat_kmeans",
            "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/20-180/mat_spectral")
    fraction = c("0-0.2","0.22-3","20-180")
    data = list(data.0_0.2$Genocenose,data.0.22_3$Genocenose,data.20_180$Genocenose)
    l = list(l1,l2,l3)
    n = list(n1,n2,n3)
    
    RI1 = matrix(0, nrow = round(n[[1]]/2), ncol = 2)
    RI2 = matrix(0, nrow = round(n[[1]]/2), ncol = 2)
    RI3 = matrix(0, nrow = round(n[[2]]/2), ncol = 2)
    RI4 = matrix(0, nrow = round(n[[2]]/2), ncol = 2)
    RI5 = matrix(0, nrow = round(n[[3]]/2), ncol = 2)
    RI6 = matrix(0, nrow = round(n[[3]]/2), ncol = 2)
    
    RI = list(RI1,RI2,RI3,RI4,RI5,RI6)
    
    for(num_RI in 1:6){
      setwd(dir=dir[num_RI])
      if(even(num_RI)){
        mat_spectral = read.table(file = matrices.list.spectral2[num_matrix2])
        mat = as.matrix(mat_spectral) 
      }else{
        mat_kmeans = read.table(file = matrices.list.kmeans2[num_matrix2])
        mat = as.matrix(mat_kmeans) 
      }
      for(k in 1:round(n[[ceiling(num_RI/2)]]/2)){
        RI[[num_RI]][k,1] = k
        RI[[num_RI]][k,2] = rand.index(mat[14*(k-1)+l[[ceiling(num_RI/2)]]-1,3:(n[[ceiling(num_RI/2)]]+2)], data[[ceiling(num_RI/2)]])
      } 
      RI[[num_RI]] = as.data.frame(RI[[num_RI]])
    }
    
    res3[[a]] <- lapply(1:6, function(i) {
      res3[[a]] <- RI[[i]]
      colnames(res3[[a]]) <- c("Param", "RI")
      res3[[a]] %>% mutate(Type = type(i), Fraction = fraction[ceiling(i/2)], Distance = select_distances_matrix)
    }) %>%
      bind_rows() 
    a = a + 1
  }
  
  res3 = bind_rows(res3[[1]], res3[[2]], res3[[3]], res3[[4]])
  
  res <- bind_rows(res1, res2, res3)
  
  res$Distance[which(res$Distance=="jaccard_abundance")] = "JA"
  res$Distance[which(res$Distance=="ab_jaccard_abundance")] = "abJ"
  res$Distance[which(res$Distance=="braycurtis_abundance")] = "BA"
  res$Distance[which(res$Distance=="ab_ochiai_abundance")] = "abO"
  res$Distance[which(res$Distance=="ab_sorensen_abundance")] = "abS"
  res$Distance[which(res$Distance=="simka_jaccard_abundance")] = "sJA"
  res$Distance[which(res$Distance=="chord_prevalence")] = "CP"
  res$Distance[which(res$Distance=="jaccard_prevalence")] = "JP"
  res$Distance[which(res$Distance=="kulczynski_prevalence")] = "KP"
  res$Distance[which(res$Distance=="ochiai_prevalence")] = "OP"
  res$Distance[which(res$Distance=="whittaker_prevalence")] = "WP"
  res$Distance[which(res$Distance=="simka_jaccard_prevalence")] = "sJP"
  res$Distance[which(res$Distance=="sorensen_braycurtis_prevalence")] = "SBP"
  res$Distance[which(res$Distance=="ochiai_abundance")] = "OA"
  res$Distance[which(res$Distance=="sorensen_abundance")] = "SA"
  res$Distance[which(res$Distance=="chord_hellinger_prevalence")] = "CHP"
  res$Distance[which(res$Distance=="jaccard_canberra_prevalence")] = "JCP"
  
  return(res)
}

stability_study_results <- function(type = NULL, method = c("RI_plot","RI_summary"), threshold = NULL, size_fraction){
  
  #################################### Representation of the RIs in function of the size of the clustering, the distances 
  #################################### matrix and the method (classical or ordinal) 
  
  if(size_fraction=="0.8-5"){
    setwd(dir="~/metaord/0.8-5/RI_results") 
  }
  
  if(size_fraction=="5-20"){
    setwd(dir="~/metaord/5-20/RI_results") 
  }
  
  if(size_fraction=="180-2000"){
    setwd(dir="~/metaord/180-2000/RI_results") 
  }
  
  if(size_fraction=="0.22-3"){
    setwd(dir="~/metaord/0.22-3/RI_results") 
  }
  
  if(size_fraction=="0-0.2"){
    setwd(dir="~/metaord/0-0.2/RI_results") 
  }
  
  if(size_fraction=="20-180"){
    setwd(dir="~/metaord/20-180/RI_results") 
  }
  if(method == "RI_plot"){ # plot of all of the RIs
    if(type == "kmeans"){ # clustering obtained by the classical method
      res = fread("RI_kmeans_all") %>% rename(param1 = d1, param2 = d2)
      res$Diff = abs(res$param1 - res$param2)
    }
    if(type == "spectral"){ # clustering obtained by the ordinal method
      res = fread("RI_spectral_all") %>% rename(param1 = k1, param2 = k2)
      res$Diff = abs(res$param1 - res$param2)
    }
    if(is.null(threshold)){ # without threshold
      p = ggplot(res, aes(x = Diff, y = RI, color = Distance)) +
        geom_point(alpha = 0.1, size = 0.5) +
        scale_y_continuous(limits = c(0,1)) +
        facet_grid(Distance~Nb_Cluster)
    }
    if(length(threshold) == 1){ # threshold on the difference between two distinct parameters (k1 - k2 or d1 - d2)
      p = ggplot(res %>% filter(Diff <= threshold), aes(x = Diff, y = RI, color = Distance)) +
        geom_point(alpha = 0.1, size = 0.5) +
        scale_y_continuous(limits = c(0,1)) +
        facet_grid(Distance~Nb_Cluster)
    }
    if(length(threshold) == 2){ # threshold directly on k1 and k2 or d1 and d2
      p = ggplot(res %>% filter(param1 >= threshold[1] & param2 >= threshold[1] & param1 <= threshold[2] & param2 <= threshold[2]), aes(x = Diff, y = RI, color = Distance)) +
        geom_point(alpha = 0.1, size = 0.5) +
        scale_y_continuous(limits = c(0,1)) +
        facet_grid(Distance~Nb_Cluster)
    }
  }
  if(method == "RI_summary"){ # plot of a summary of the RIs for each method : median, q10 and q90 
    res1 = fread("RI_kmeans_all") %>% rename(param1 = d1, param2 = d2)
    res2 = fread("RI_spectral_all") %>% rename(param1 = k1, param2 = k2)
    res <- bind_rows(res1 %>% mutate(Type = "K-means"), 
                     res2 %>% mutate(Type = "Spectral")) %>% 
      mutate(Diff = abs(param1 - param2))
    if(is.null(threshold)){
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
    }
    if(length(threshold) == 1){
      res1 <- res %>% filter(Diff <= threshold) %>% group_by(Nb_Cluster, Distance, Diff, Type) %>% 
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
    }
    if(length(threshold) == 2){
      res1 <- res %>% filter(param1 >= threshold[1] & param2 >= threshold[1] & param1 <= threshold[2] & param2 <= threshold[2]) %>% group_by(Nb_Cluster, Distance, Diff, Type) %>% 
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
    }
  }
  return(p)
}

spectral.clustering.new <- function(A, normalised = TRUE, score = FALSE, K = 2, adj = FALSE){
  
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
  U = scale(U, center = F)
  temp = unique(U, margin = 2)
  if(dim(temp)[1] < K){stop('FAIL!')}
  k.means = kmeans(U, centers = K, nstart = 1000, iter.max = 20)
  labels[noniso.seq] = k.means$cluster
  
  return(labels)
}

identify_samples <- function(size_fraction){
  
  if(size_fraction == "0-0.2"){
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
    
    K_means_initial = K_means
    all_permutations = t(array(unlist(permn(1:l)), dim = c(l, fact(l))))
    score = 0
    score_ex = 0
    a = 0
    K_means_ex = K_means
    
    for(j in 1:fact(l)){
      for(perm_l in 1:l){
        K_means[which(K_means_ex==perm_l)] = all_permutations[j,perm_l]
      }
      score_ex = sum(diag(table(K_means,h_clust)))
      if(score_ex >= score){
        score = score_ex
        a = j
      }
      K_means = K_means_ex
    }
    
    for(perm_l in 1:l){
      K_means[which(K_means_ex==perm_l)] = all_permutations[a,perm_l]
    }
  }
  
  #############################################################################################
  
  if(size_fraction == "0.22-3"){
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
    
    K_means_initial = K_means
    all_permutations = t(array(unlist(permn(1:l)), dim = c(l, fact(l))))
    score = 0
    score_ex = 0
    a = 0
    K_means_ex = K_means
    
    for(j in 1:fact(l)){
      for(perm_l in 1:l){
        K_means[which(K_means_ex==perm_l)] = all_permutations[j,perm_l]
      }
      score_ex = sum(diag(table(K_means,h_clust)))
      if(score_ex >= score){
        score = score_ex
        a = j
      }
      K_means = K_means_ex
    }
    
    for(perm_l in 1:l){
      K_means[which(K_means_ex==perm_l)] = all_permutations[a,perm_l]
    }
  }
  
  #############################################################################################
  
  if(size_fraction == "0.8-5"){
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
    
    K_means_initial = K_means
    all_permutations = t(array(unlist(permn(1:l)), dim = c(l, fact(l))))
    score = 0
    score_ex = 0
    a = 0
    K_means_ex = K_means
    
    for(j in 1:fact(l)){
      for(perm_l in 1:l){
        K_means[which(K_means_ex==perm_l)] = all_permutations[j,perm_l]
      }
      score_ex = sum(diag(table(K_means,h_clust)))
      if(score_ex >= score){
        score = score_ex
        a = j
      }
      K_means = K_means_ex
    }
    
    for(perm_l in 1:l){
      K_means[which(K_means_ex==perm_l)] = all_permutations[a,perm_l]
    }
  }
  
  
  #############################################################################################
  
  if(size_fraction == "5-20"){
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
    
    K_means_initial = K_means
    all_permutations = t(array(unlist(permn(1:l)), dim = c(l, fact(l))))
    score = 0
    score_ex = 0
    a = 0
    K_means_ex = K_means
    
    for(j in 1:fact(l)){
      for(perm_l in 1:l){
        K_means[which(K_means_ex==perm_l)] = all_permutations[j,perm_l]
      }
      score_ex = sum(diag(table(K_means,h_clust)))
      if(score_ex >= score){
        score = score_ex
        a = j
      }
      K_means = K_means_ex
    }
    
    for(perm_l in 1:l){
      K_means[which(K_means_ex==perm_l)] = all_permutations[a,perm_l]
    }
  }
  
  #############################################################################################
  
  if(size_fraction == "20-180"){
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
    
    K_means_initial = K_means
    all_permutations = t(array(unlist(permn(1:l)), dim = c(l, fact(l))))
    score = 0
    score_ex = 0
    a = 0
    K_means_ex = K_means
    
    for(j in 1:fact(l)){
      for(perm_l in 1:l){
        K_means[which(K_means_ex==perm_l)] = all_permutations[j,perm_l]
      }
      score_ex = sum(diag(table(K_means,h_clust)))
      if(score_ex >= score){
        score = score_ex
        a = j
      }
      K_means = K_means_ex
    }
    
    for(perm_l in 1:l){
      K_means[which(K_means_ex==perm_l)] = all_permutations[a,perm_l]
    }
  }
  
  #############################################################################################
  
  if(size_fraction == "180-2000"){
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
    
    K_means_initial = K_means
    all_permutations = t(array(unlist(permn(1:l)), dim = c(l, fact(l))))
    score = 0
    score_ex = 0
    a = 0
    K_means_ex = K_means
    
    for(j in 1:fact(l)){
      for(perm_l in 1:l){
        K_means[which(K_means_ex==perm_l)] = all_permutations[j,perm_l]
      }
      score_ex = sum(diag(table(K_means,h_clust)))
      if(score_ex >= score){
        score = score_ex
        a = j
      }
      K_means = K_means_ex
    }
    
    for(perm_l in 1:l){
      K_means[which(K_means_ex==perm_l)] = all_permutations[a,perm_l]
    }
  }
  
  return(list(Kmeans_clustering=K_means,Genoscope_clustering=h_clust,Kmeans_initial=K_means_initial,score=score))
}

eigenvalues_study <- function(size_fraction){
  
  matrices.list1 = c("jaccard_abundance", "ab_jaccard_abundance", "braycurtis_abundance", 
                     "ab_ochiai_abundance", "ab_sorensen_abundance", "simka_jaccard_abundance", 
                     "chord_prevalence", "jaccard_prevalence", "kulczynski_prevalence", 
                     "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence",
                     "sorensen_braycurtis_prevalence")
  matrices.list2 = c("jaccard_abundance", "ochiai_abundance", "sorensen_abundance", 
                     "simka_jaccard_abundance", "chord_hellinger_prevalence", 
                     "jaccard_canberra_prevalence", "kulczynski_prevalence", 
                     "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence",
                     "sorensen_braycurtis_prevalence")
  
  if(size_fraction == "0-0.2"){
    delete = which(is.na(data.0_0.2$Genocenose))
    data.0_0.2 = data.0_0.2[-delete,]
    import_data(size_fraction, samples = rownames(data.0_0.2))
    data.0_0.2 = data.0_0.2[metagenomic_sample,]
    n = dim(jaccard_abundance)[1]
    l = 8
    mds_eig = matrix(0,length(matrices.list2),n+1)
    
    a = 1
    for(mat in matrices.list2){
      mds_eig[a,1] = a
      mds_eig[a,2:(n+1)] = cmdscale(get(mat), k = 3, eig = T)$eig
      a = a + 1
    }
  }
  
  #############################################################################################
  
  if(size_fraction == "0.22-3"){
    delete = which(is.na(data.0.22_3$Genocenose))
    data.0.22_3 = data.0.22_3[-delete,]
    import_data(size_fraction, samples = rownames(data.0.22_3))
    data.0.22_3 = data.0.22_3[metagenomic_sample,]
    n = dim(jaccard_abundance)[1]
    l = 8
    mds_eig = matrix(0,length(matrices.list2),n+1)
    
    a = 1
    for(mat in matrices.list2){
      mds_eig[a,1] = a
      mds_eig[a,2:(n+1)] = cmdscale(get(mat), k = 3, eig = T)$eig
      a = a + 1
    }
  }
  
  #############################################################################################
  
  if(size_fraction == "0.8-5"){
    delete = which(is.na(data.0.8_5$Genocenose))
    data.0.8_5 = data.0.8_5[-delete,]
    import_data(size_fraction, samples = rownames(data.0.8_5))
    data.0.8_5 = data.0.8_5[metagenomic_sample,]
    n = dim(jaccard_abundance)[1]
    l = 11
    mds_eig = matrix(0,length(matrices.list1),n+1)
    
    a = 1
    for(mat in matrices.list1){
      mds_eig[a,1] = a
      mds_eig[a,2:(n+1)] = cmdscale(get(mat), k = 3, eig = T)$eig
      a = a + 1
    }
  }
  
  
  #############################################################################################
  
  if(size_fraction == "5-20"){
    delete = which(is.na(data.5_20$Genocenose))
    data.5_20 = data.5_20[-delete,]
    import_data(size_fraction, samples = rownames(data.5_20))
    data.5_20 = data.5_20[metagenomic_sample,]
    n = dim(jaccard_abundance)[1]
    l = 6
    mds_eig = matrix(0,length(matrices.list1),n+1)
    
    a = 1
    for(mat in matrices.list1){
      mds_eig[a,1] = a
      mds_eig[a,2:(n+1)] = cmdscale(get(mat), k = 3, eig = T)$eig
      a = a + 1
    }
  }
  
  #############################################################################################
  
  if(size_fraction == "20-180"){
    delete = which(is.na(data.20_180$Genocenose))
    data.20_180 = data.20_180[-delete,]
    import_data(size_fraction, samples = rownames(data.20_180))
    data.20_180 = data.20_180[metagenomic_sample,]
    n = dim(jaccard_abundance)[1]
    l = 6
    mds_eig = matrix(0,length(matrices.list2),n+1)
    
    a = 1
    for(mat in matrices.list2){
      mds_eig[a,1] = a
      mds_eig[a,2:(n+1)] = cmdscale(get(mat), k = 3, eig = T)$eig
      a = a + 1
    }
  }
    
    if(size_fraction == "180-2000"){
      import_data(size_fraction, samples = rownames(data.180_2000))
      data.180_2000 = data.180_2000[metagenomic_sample,]
      n = dim(jaccard_abundance)[1]
      l = 8
      mds_eig = matrix(0,length(matrices.list1),n+1)
      
      a = 1
      for(mat in matrices.list1){
        mds_eig[a,1] = a
        mds_eig[a,2:(n+1)] = cmdscale(get(mat), k = 3, eig = T)$eig
        a = a + 1
      }
    }
  
  mat_names1 = c("abJ", "abO", "abS", "BA", "CP", "JA", "JP", "KP", "OP", "SBP", "sJA", "sJP", "WP")
  mat_names2 = c("CHP", "JA", "JCP", "KP", "OA", "OP", "SA", "SBP", "sJA", "sJP", "WP")
  
  if(size_fraction=="0.8-5" | size_fraction=="5-20" | size_fraction=="180-2000"){
    par(mfrow=c(4,4))
    for(a in 1:length(matrices.list1)){
      plot(1:n,mds_eig[a,2:(n+1)], xlab = "Dim", ylab = "v.p", main = mat_names1[a], col = "blue", pch = 20) 
    }
  }else{
    par(mfrow=c(3,4))
    for(a in 1:length(matrices.list2)){
      plot(1:n,mds_eig[a,2:(n+1)], xlab = "Dim", ylab = "v.p", main = mat_names2[a], col = "blue", pch = 20) 
    }
  }
}

# identify_samples_ex <- function(size_fraction){
#   
#   if(size_fraction == "0-0.2"){
#     dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0-0.2/mat_kmeans")
#     setwd(dir=dir)
#     l = max(data.0_0.2$Genocenose)
#     index = which(comparison$Distance=="SBP" & comparison$Fraction==size_fraction & comparison$Type=="Kmeans")
#     new_RI = comparison[index,]$RI
#     k = which.max(new_RI)
#     RI_max = new_RI[k]
#     mat = read.table("Mat11_kmeans")
#     
#     h_clust = data.0_0.2$Genocenose
#     n = dim(data.0_0.2)[1]
#     K_means = rep(NA,n)
#     
#     for(i in 3:(n+2)){
#       K_means[i-2] = mat[which(mat[,1]==k & mat[,2]==l),i]
#     }
#     
#     K_means_initial = K_means
#     all_permutations = t(array(unlist(permn(1:l)), dim = c(l, fact(l))))
#     score = rep(0,fact(l))
#     score_ex = 0
#     a = 0
#     K_means_ex = K_means
#     
#     for(j in 1:fact(l)){
#       for(perm_l in 1:l){
#         K_means[which(K_means_ex==perm_l)] = all_permutations[j,perm_l] 
#       }
#       score_ex = sum(diag(table(K_means,h_clust)))
#       score[j] = score_ex
#       a = j
#       K_means = K_means_ex 
#     }    
#     
#     for(perm_l in 1:l){
#       K_means[which(K_means_ex==perm_l)] = all_permutations[a,perm_l] 
#     }
#   } 
#   
#   #############################################################################################
#   
#   if(size_fraction == "0.22-3"){
#     dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.22-3/mat_kmeans")
#     setwd(dir=dir)
#     l = max(data.0.22_3$Genocenose)
#     index = which(comparison$Distance=="SBP" & comparison$Fraction==size_fraction & comparison$Type=="Kmeans")
#     new_RI = comparison[index,]$RI
#     k = which.max(new_RI)
#     RI_max = new_RI[k]
#     which(matrices.list=="sorensen_braycurtis_prevalence")
#     mat = read.table("Mat11_kmeans")
#     
#     h_clust = data.0.22_3$Genocenose
#     n = dim(data.0.22_3)[1]
#     K_means = rep(NA,n)
#     
#     for(i in 3:(n+2)){
#       K_means[i-2] = mat[which(mat[,1]==k & mat[,2]==l),i]
#     } 
#     
#     K_means_initial = K_means
#     all_permutations = t(array(unlist(permn(1:l)), dim = c(l, fact(l))))
#     score = rep(0,fact(l))
#     score_ex = 0
#     a = 0
#     K_means_ex = K_means
#     
#     for(j in 1:fact(l)){
#       for(perm_l in 1:l){
#         K_means[which(K_means_ex==perm_l)] = all_permutations[j,perm_l] 
#       }
#       score_ex = sum(diag(table(K_means,h_clust)))
#       score[j] = score_ex
#       a = j
#       K_means = K_means_ex 
#     }    
#     
#     for(perm_l in 1:l){
#       K_means[which(K_means_ex==perm_l)] = all_permutations[a,perm_l] 
#     } 
#   }
#   
#   #############################################################################################
#   
#   if(size_fraction == "0.8-5"){
#     dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.8-5/mat_kmeans")
#     setwd(dir=dir)
#     l = max(data.0.8_5$Genocenose)
#     index = which(comparison$Distance=="SBP" & comparison$Fraction==size_fraction & comparison$Type=="Kmeans")
#     new_RI = comparison[index,]$RI
#     k = which.max(new_RI)
#     RI_max = new_RI[k]
#     which(matrices.list=="sorensen_braycurtis_prevalence")
#     mat = read.table("Mat13_kmeans")
#     
#     h_clust = data.0.8_5$Genocenose
#     n = dim(data.0.8_5)[1]
#     K_means = rep(NA,n)
#     
#     for(i in 3:(n+2)){
#       K_means[i-2] = mat[which(mat[,1]==k & mat[,2]==l),i]
#     }    
#     
#     K_means_initial = K_means
#     all_permutations = t(array(unlist(permn(1:l)), dim = c(l, fact(l))))
#     score = rep(0,fact(l))
#     score_ex = 0
#     a = 0
#     K_means_ex = K_means
#     
#     for(j in 1:fact(l)){
#       for(perm_l in 1:l){
#         K_means[which(K_means_ex==perm_l)] = all_permutations[j,perm_l] 
#       }
#       score_ex = sum(diag(table(K_means,h_clust)))
#       score[j] = score_ex
#       a = j
#       K_means = K_means_ex 
#       print(j)
#     }    
#     
#     for(perm_l in 1:l){
#       K_means[which(K_means_ex==perm_l)] = all_permutations[a,perm_l] 
#     }   
#   }
#   
#   
#   #############################################################################################
#   
#   if(size_fraction == "5-20"){
#     dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/5-20/mat_kmeans")
#     setwd(dir=dir)
#     l = max(data.5_20$Genocenose)
#     index = which(comparison$Distance=="SBP" & comparison$Fraction==size_fraction & comparison$Type=="Kmeans")
#     new_RI = comparison[index,]$RI
#     k = which.max(new_RI)
#     RI_max = new_RI[k]
#     which(matrices.list=="sorensen_braycurtis_prevalence")
#     mat = read.table("Mat13_kmeans")
#     
#     h_clust = data.5_20$Genocenose
#     n = dim(data.5_20)[1]
#     K_means = rep(NA,n)
#     
#     for(i in 3:(n+2)){
#       K_means[i-2] = mat[which(mat[,1]==k & mat[,2]==l),i]
#     } 
#     
#     K_means_initial = K_means
#     all_permutations = t(array(unlist(permn(1:l)), dim = c(l, fact(l))))
#     score = rep(0,fact(l))
#     score_ex = 0
#     a = 0
#     K_means_ex = K_means
#     
#     for(j in 1:fact(l)){
#       for(perm_l in 1:l){
#         K_means[which(K_means_ex==perm_l)] = all_permutations[j,perm_l] 
#       }
#       score_ex = sum(diag(table(K_means,h_clust)))
#       score[j] = score_ex
#       a = j
#       K_means = K_means_ex 
#     }    
#     
#     for(perm_l in 1:l){
#       K_means[which(K_means_ex==perm_l)] = all_permutations[a,perm_l] 
#     }  
#   }
#   
#   #############################################################################################
#   
#   if(size_fraction == "20-180"){
#     dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/20-180/mat_kmeans")
#     setwd(dir=dir)
#     l = max(data.20_180$Genocenose)
#     index = which(comparison$Distance=="SBP" & comparison$Fraction==size_fraction & comparison$Type=="Kmeans")
#     new_RI = comparison[index,]$RI
#     k = which.max(new_RI)
#     RI_max = new_RI[k]
#     which(matrices.list=="sorensen_braycurtis_prevalence")
#     mat = read.table("Mat11_kmeans")
#     
#     h_clust = data.20_180$Genocenose
#     n = dim(data.20_180)[1]
#     K_means = rep(NA,n)
#     
#     for(i in 3:(n+2)){
#       K_means[i-2] = mat[which(mat[,1]==k & mat[,2]==l),i]
#     } 
#     
#     K_means_initial = K_means
#     all_permutations = t(array(unlist(permn(1:l)), dim = c(l, fact(l))))
#     score = rep(0,fact(l))
#     score_ex = 0
#     a = 0
#     K_means_ex = K_means
#     
#     for(j in 1:fact(l)){
#       for(perm_l in 1:l){
#         K_means[which(K_means_ex==perm_l)] = all_permutations[j,perm_l] 
#       }
#       score_ex = sum(diag(table(K_means,h_clust)))
#       score[j] = score_ex
#       a = j
#       K_means = K_means_ex 
#     }    
#     
#     for(perm_l in 1:l){
#       K_means[which(K_means_ex==perm_l)] = all_permutations[a,perm_l] 
#     } 
#   }
#   
#   #############################################################################################
#   
#   if(size_fraction == "180-2000"){
#     dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/180-2000/mat_kmeans")
#     setwd(dir=dir)
#     l = max(data.180_2000$Genocenose)
#     index = which(comparison$Distance=="SBP" & comparison$Fraction==size_fraction & comparison$Type=="Kmeans")
#     new_RI = comparison[index,]$RI
#     k = which.max(new_RI)
#     RI_max = new_RI[k]
#     which(matrices.list=="sorensen_braycurtis_prevalence")
#     mat = read.table("Mat13_kmeans")
#     
#     h_clust = data.180_2000$Genocenose
#     n = dim(data.180_2000)[1]
#     K_means = rep(NA,n)
#     
#     for(i in 3:(n+2)){
#       K_means[i-2] = mat[which(mat[,1]==k & mat[,2]==l),i]
#     }
#     
#     K_means_initial = K_means
#     all_permutations = t(array(unlist(permn(1:l)), dim = c(l, fact(l))))
#     score = rep(0,fact(l))
#     score_ex = 0
#     a = 0
#     K_means_ex = K_means
#     
#     for(j in 1:fact(l)){
#       for(perm_l in 1:l){
#         K_means[which(K_means_ex==perm_l)] = all_permutations[j,perm_l] 
#       }
#       score_ex = sum(diag(table(K_means,h_clust)))
#       score[j] = score_ex
#       a = j
#       K_means = K_means_ex 
#     }    
#     
#     for(perm_l in 1:l){
#       K_means[which(K_means_ex==perm_l)] = all_permutations[a,perm_l] 
#     } 
#   }
#   
#   return(list(all_permutations=all_permutations,score=score))
# }



# compare_clusterings_ex <- function(){
# 
#   ########################### Function comparing distinct clusterings by the calculation of the Rand index
# 
#   library(gtools)
# 
#   matrices.list1 = c("jaccard_abundance", "ab_jaccard_abundance", "braycurtis_abundance",
#                      "ab_ochiai_abundance", "ab_sorensen_abundance", "simka_jaccard_abundance",
#                      "chord_prevalence", "jaccard_prevalence", "kulczynski_prevalence",
#                      "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence",
#                      "sorensen_braycurtis_prevalence")
#   matrices.list2 = c("jaccard_abundance", "ochiai_abundance", "sorensen_abundance",
#                      "simka_jaccard_abundance", "chord_hellinger_prevalence",
#                      "jaccard_canberra_prevalence", "kulczynski_prevalence",
#                      "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence",
#                      "sorensen_braycurtis_prevalence")
#   matrices.list.spectral1 = c("Mat1_spectral", "Mat2_spectral", "Mat3_spectral", "Mat4_spectral",
#                               "Mat5_spectral", "Mat6_spectral", "Mat7_spectral", "Mat8_spectral",
#                               "Mat9_spectral", "Mat10_spectral", "Mat11_spectral", "Mat12_spectral",
#                               "Mat13_spectral")
#   matrices.list.kmeans1 = c("Mat1_kmeans", "Mat2_kmeans", "Mat3_kmeans", "Mat4_kmeans",
#                             "Mat5_kmeans", "Mat6_kmeans", "Mat7_kmeans", "Mat8_kmeans",
#                             "Mat9_kmeans", "Mat10_kmeans", "Mat11_kmeans", "Mat12_kmeans",
#                             "Mat13_kmeans")
# 
#   matrices.list.spectral2  = c("Mat1_spectral", "Mat2_spectral", "Mat3_spectral", "Mat4_spectral",
#                                "Mat5_spectral", "Mat6_spectral", "Mat7_spectral", "Mat8_spectral",
#                                "Mat9_spectral", "Mat10_spectral", "Mat11_spectral")
#   matrices.list.kmeans2    = c("Mat1_kmeans", "Mat2_kmeans", "Mat3_kmeans", "Mat4_kmeans",
#                                "Mat5_kmeans", "Mat6_kmeans", "Mat7_kmeans", "Mat8_kmeans",
#                                "Mat9_kmeans", "Mat10_kmeans", "Mat11_kmeans")
# 
#   size_fraction = "0.8-5"
#   delete = which(is.na(data.0.8_5$Genocenose))
#   data.0.8_5 = data.0.8_5[-delete,]
#   import_data(size_fraction, samples = rownames(data.0.8_5))
#   data.0.8_5 = data.0.8_5[metagenomic_sample,]
# 
#   size_fraction = "0-0.2"
#   delete = which(is.na(data.0_0.2$Genocenose))
#   data.0_0.2 = data.0_0.2[-delete,]
#   import_data(size_fraction, samples = rownames(data.0_0.2))
#   data.0_0.2 = data.0_0.2[metagenomic_sample,]
# 
#   size_fraction = "0.22-3"
#   delete = which(is.na(data.0.22_3$Genocenose))
#   data.0.22_3 = data.0.22_3[-delete,]
#   import_data(size_fraction, samples = rownames(data.0.22_3))
#   data.0.22_3 = data.0.22_3[metagenomic_sample,]
# 
#   size_fraction = "5-20"
#   delete = which(is.na(data.5_20$Genocenose))
#   data.5_20 = data.5_20[-delete,]
#   import_data(size_fraction, samples = rownames(data.5_20))
#   data.5_20 = data.5_20[metagenomic_sample,]
# 
#   size_fraction = "20-180"
#   delete = which(is.na(data.20_180$Genocenose))
#   data.20_180 = data.20_180[-delete,]
#   import_data(size_fraction, samples = rownames(data.20_180))
#   data.20_180 = data.20_180[metagenomic_sample,]
# 
#   size_fraction = "180-2000"
#   # delete = which(is.na(data.180_2000$Genocenose))
#   # data.180_2000 = data.180_2000[-delete,]
#   import_data(size_fraction, samples = rownames(data.180_2000))
#   data.180_2000 = data.180_2000[metagenomic_sample,]
# 
#   res_1 = data.frame()
#   res_2 = data.frame()
#   res_3 = data.frame()
#   res_4 = data.frame()
#   res_5 = data.frame()
#   res_6 = data.frame()
#   res_7 = data.frame()
# 
#   res1 = list(res_1,res_2,res_3,res_4,res_5,res_6,res_7)
#   res2 = list(res_1,res_2,res_3,res_4,res_5,res_6)
#   res3 = list(res_1,res_2,res_3,res_4)
# 
#   res = list(res1,res2,res3)
# 
#   a = 1
# 
#   ################################### Compare the clusterings obtained by Genoscope and those obtained
#   ################################### by the two methods (classical and ordinal) for distances matrix
#   ################################### which are present in all of the 6 sizes fraction
# 
#   for(select_distances_matrix in c("jaccard_abundance", "simka_jaccard_abundance",
#                                    "kulczynski_prevalence", "ochiai_prevalence",
#                                    "whittaker_prevalence", "simka_jaccard_prevalence",
#                                    "sorensen_braycurtis_prevalence")){
#     num_matrix1 = which(matrices.list1 == select_distances_matrix)
#     num_matrix2 = which(matrices.list2 == select_distances_matrix)
# 
#     l1 = max(data.0.8_5$Genocenose)
#     l2 = max(data.0_0.2$Genocenose)
#     l3 = max(data.0.22_3$Genocenose)
#     l4 = max(data.5_20$Genocenose)
#     l5 = max(data.20_180$Genocenose)
#     l6 = max(data.180_2000$Genocenose)
# 
#     n1 = dim(data.0.8_5)[1]
#     n3 = dim(data.0_0.2)[1]
#     n5 = dim(data.0.22_3)[1]
#     n7 = dim(data.5_20)[1]
#     n9 = dim(data.20_180)[1]
#     n11 = dim(data.180_2000)[1]
#     n2 = dim(data.0.8_5)[1]
#     n4 = dim(data.0_0.2)[1]
#     n6 = dim(data.0.22_3)[1]
#     n8 = dim(data.5_20)[1]
#     n10 = dim(data.20_180)[1]
#     n12 = dim(data.180_2000)[1]
# 
#     dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.8-5/mat_kmeans",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.8-5/mat_spectral",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0-0.2/mat_kmeans",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0-0.2/mat_spectral",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.22-3/mat_kmeans",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.22-3/mat_spectral",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/5-20/mat_kmeans",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/5-20/mat_spectral",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/20-180/mat_kmeans",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/20-180/mat_spectral",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/180-2000/mat_kmeans",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/180-2000/mat_spectral")
#     type = rep(c("Kmeans","Spectral"),6)
#     fraction = c("0.8-5","0.8-5","0-0.2","0-0.2","0.22-3","0.22-3","5-20","5-20","20-180","20-180","180-2000","180-2000")
#     data = list(data.0.8_5$Genocenose,data.0_0.2$Genocenose,data.0.22_3$Genocenose,data.5_20$Genocenose,
#                 data.20_180$Genocenose,data.180_2000$Genocenose)
#     l = list(l1,l2,l3,l4,l5,l6)
#     n = list(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12)
# 
#     RI1 = matrix(0, nrow = round(n[[1]]/2), ncol = 2)
#     RI2 = matrix(0, nrow = round(n[[2]]/2), ncol = 2)
#     RI3 = matrix(0, nrow = round(n[[3]]/2), ncol = 2)
#     RI4 = matrix(0, nrow = round(n[[4]]/2), ncol = 2)
#     RI5 = matrix(0, nrow = round(n[[5]]/2), ncol = 2)
#     RI6 = matrix(0, nrow = round(n[[6]]/2), ncol = 2)
#     RI7 = matrix(0, nrow = round(n[[7]]/2), ncol = 2)
#     RI8 = matrix(0, nrow = round(n[[8]]/2), ncol = 2)
#     RI9 = matrix(0, nrow = round(n[[9]]/2), ncol = 2)
#     RI10 = matrix(0, nrow = round(n[[10]]/2), ncol = 2)
#     RI11 = matrix(0, nrow = round(n[[11]]/2), ncol = 2)
#     RI12 = matrix(0, nrow = round(n[[12]]/2), ncol = 2)
# 
#     RI = list(RI1,RI2,RI3,RI4,RI5,RI6,RI7,RI8,RI9,RI10,RI11,RI12)
# 
#     for(num_RI in 1:12){
#       setwd(dir=dir[num_RI])
#       if(even(num_RI) & num_RI == 2 | num_RI == 8 |num_RI == 12){
#         mat_spectral = read.table(file = matrices.list.spectral1[num_matrix1])
#         mat = as.matrix(mat_spectral)
#       }
#       if(even(num_RI) & num_RI == 4 | num_RI == 6 |num_RI == 10){
#         mat_spectral = read.table(file = matrices.list.spectral2[num_matrix2])
#         mat = as.matrix(mat_spectral)
#       }
#       if(odd(num_RI) & num_RI == 1 | num_RI == 7 |num_RI == 11){
#         mat_kmeans = read.table(file = matrices.list.kmeans1[num_matrix1])
#         mat = as.matrix(mat_kmeans)
#       }
#       if(odd(num_RI) & num_RI == 3 | num_RI == 5 |num_RI == 9){
#         mat_kmeans = read.table(file = matrices.list.kmeans2[num_matrix2])
#         mat = as.matrix(mat_kmeans)
#       }
#       for(k in 1:round(n[[num_RI]]/2)){
#         RI[[num_RI]][k,1] = k
#         RI[[num_RI]][k,2] = rand.index(mat[14*(k-1)+l[[ceiling(num_RI/2)]]-1,3:(n[[num_RI]]+2)], data[[ceiling(num_RI/2)]])
#       }
#       RI[[num_RI]] = as.data.frame(RI[[num_RI]])
#     }
# 
#     res1[[a]] <- lapply(1:12, function(i) {
#       res1[[a]] <- RI[[i]]
#       colnames(res1[[a]]) <- c("Param", "RI")
#       res1[[a]] %>% mutate(Type = type[i], Fraction = fraction[i], Distance = select_distances_matrix)
#     }) %>%
#       bind_rows()
#     a = a + 1
#   }
# 
#   res1 = bind_rows(res1[[1]], res1[[2]], res1[[3]], res1[[4]], res1[[5]], res1[[6]], res1[[7]])
#   a = 1
# 
#   ################################### Compare the clusterings obtained by Genoscope and those obtained
#   ################################### by the two methods (classical and ordinal) for distances matrix
#   ################################### which are only present in 3 of the 6 sizes fraction
# 
#   for(select_distances_matrix in c("ab_jaccard_abundance", "braycurtis_abundance",
#                                    "ab_ochiai_abundance", "ab_sorensen_abundance",
#                                    "chord_prevalence", "jaccard_prevalence")){
#     num_matrix1 = which(matrices.list1 == select_distances_matrix)
# 
#     l1 = max(data.0.8_5$Genocenose)
#     l2 = max(data.5_20$Genocenose)
#     l3 = max(data.180_2000$Genocenose)
# 
#     n1 = dim(data.0.8_5)[1]
#     n3 = dim(data.5_20)[1]
#     n5 = dim(data.180_2000)[1]
#     n2 = dim(data.0.8_5)[1]
#     n4 = dim(data.5_20)[1]
#     n6 = dim(data.180_2000)[1]
# 
#     dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.8-5/mat_kmeans",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.8-5/mat_spectral",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/5-20/mat_kmeans",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/5-20/mat_spectral",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/180-2000/mat_kmeans",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/180-2000/mat_spectral")
#     type = rep(c("Kmeans","Spectral"),3)
#     fraction = c("0.8-5","0.8-5","5-20","5-20","180-2000","180-2000")
#     data = list(data.0.8_5$Genocenose,data.5_20$Genocenose,data.180_2000$Genocenose)
#     l = list(l1,l2,l3)
#     n = list(n1,n2,n3,n4,n5,n6)
# 
#     RI1 = matrix(0, nrow = round(n[[1]]/2), ncol = 2)
#     RI2 = matrix(0, nrow = round(n[[2]]/2), ncol = 2)
#     RI3 = matrix(0, nrow = round(n[[3]]/2), ncol = 2)
#     RI4 = matrix(0, nrow = round(n[[4]]/2), ncol = 2)
#     RI5 = matrix(0, nrow = round(n[[5]]/2), ncol = 2)
#     RI6 = matrix(0, nrow = round(n[[6]]/2), ncol = 2)
# 
#     RI = list(RI1,RI2,RI3,RI4,RI5,RI6)
# 
#     for(num_RI in 1:6){
#       setwd(dir=dir[num_RI])
#       if(even(num_RI)){
#         mat_spectral = read.table(file = matrices.list.spectral1[num_matrix1])
#         mat = as.matrix(mat_spectral)
#       }else{
#         mat_kmeans = read.table(file = matrices.list.kmeans1[num_matrix1])
#         mat = as.matrix(mat_kmeans)
#       }
#       for(k in 1:round(n[[num_RI]]/2)){
#         RI[[num_RI]][k,1] = k
#         RI[[num_RI]][k,2] = rand.index(mat[14*(k-1)+l[[ceiling(num_RI/2)]]-1,3:(n[[num_RI]]+2)], data[[ceiling(num_RI/2)]])
#       }
#       RI[[num_RI]] = as.data.frame(RI[[num_RI]])
#     }
# 
#     res2[[a]] <- lapply(1:6, function(i) {
#       res2[[a]] <- RI[[i]]
#       colnames(res2[[a]]) <- c("Param", "RI")
#       res2[[a]] %>% mutate(Type = type[i], Fraction = fraction[i], Distance = select_distances_matrix)
#     }) %>%
#       bind_rows()
#     a = a + 1
#   }
# 
#   res2 = bind_rows(res2[[1]], res2[[2]], res2[[3]], res2[[4]], res2[[5]], res2[[6]])
#   a = 1
# 
#   ################################### Compare the clusterings obtained by Genoscope and those obtained
#   ################################### by the two methods (classical and ordinal) for distances matrix
#   ################################### which are only present in the 3 remaining sizes fraction
# 
#   for(select_distances_matrix in c("ochiai_abundance", "sorensen_abundance", "chord_hellinger_prevalence",
#                                    "jaccard_canberra_prevalence")){
#     num_matrix2 = which(matrices.list2 == select_distances_matrix)
# 
#     l1 = max(data.0_0.2$Genocenose)
#     l2 = max(data.0.22_3$Genocenose)
#     l3 = max(data.20_180$Genocenose)
# 
#     n1 = dim(data.0_0.2)[1]
#     n3 = dim(data.0.22_3)[1]
#     n5 = dim(data.20_180)[1]
#     n2 = dim(data.0_0.2)[1]
#     n4 = dim(data.0.22_3)[1]
#     n6 = dim(data.20_180)[1]
# 
#     dir = c("~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0-0.2/mat_kmeans",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0-0.2/mat_spectral",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.22-3/mat_kmeans",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/0.22-3/mat_spectral",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/20-180/mat_kmeans",
#             "~/Bureau/CDDrnarci/stability_study/cluster_metaord_subsample_stations/20-180/mat_spectral")
#     type = rep(c("Kmeans","Spectral"),3)
#     fraction = c("0-0.2","0-0.2","0.22-3","0.22-3","20-180","20-180")
#     data = list(data.0_0.2$Genocenose,data.0.22_3$Genocenose,data.20_180$Genocenose)
#     l = list(l1,l2,l3)
#     n = list(n1,n2,n3,n4,n5,n6)
# 
#     RI1 = matrix(0, nrow = round(n[[1]]/2), ncol = 2)
#     RI2 = matrix(0, nrow = round(n[[2]]/2), ncol = 2)
#     RI3 = matrix(0, nrow = round(n[[3]]/2), ncol = 2)
#     RI4 = matrix(0, nrow = round(n[[4]]/2), ncol = 2)
#     RI5 = matrix(0, nrow = round(n[[5]]/2), ncol = 2)
#     RI6 = matrix(0, nrow = round(n[[6]]/2), ncol = 2)
# 
#     RI = list(RI1,RI2,RI3,RI4,RI5,RI6)
# 
#     for(num_RI in 1:6){
#       setwd(dir=dir[num_RI])
#       if(even(num_RI)){
#         mat_spectral = read.table(file = matrices.list.spectral2[num_matrix2])
#         mat = as.matrix(mat_spectral)
#       }else{
#         mat_kmeans = read.table(file = matrices.list.kmeans2[num_matrix2])
#         mat = as.matrix(mat_kmeans)
#       }
#       for(k in 1:round(n[[num_RI]]/2)){
#         RI[[num_RI]][k,1] = k
#         RI[[num_RI]][k,2] = rand.index(mat[14*(k-1)+l[[ceiling(num_RI/2)]]-1,3:(n[[num_RI]]+2)], data[[ceiling(num_RI/2)]])
#       }
#       RI[[num_RI]] = as.data.frame(RI[[num_RI]])
#     }
# 
#     res3[[a]] <- lapply(1:6, function(i) {
#       res3[[a]] <- RI[[i]]
#       colnames(res3[[a]]) <- c("Param", "RI")
#       res3[[a]] %>% mutate(Type = type[i], Fraction = fraction[i], Distance = select_distances_matrix)
#     }) %>%
#       bind_rows()
#     a = a + 1
#   }
# 
#   res3 = bind_rows(res3[[1]], res3[[2]], res3[[3]], res3[[4]])
# 
#   res <- bind_rows(res1, res2, res3)
# 
#   return(res)
# }

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

# Reciprocal <- function(x){
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
# laplacian <- function(A, normalised = FALSE){
# 
#     n = dim(A)[1]
#     temp = apply(abs(A), 2, sum)
#     D = diag(temp, nrow = n)
# 
#     temp1 = Reciprocal(sqrt(temp))
#     half.D = diag(temp1, nrow = n)
#     if(normalised == TRUE) 	return(half.D %*% (D - A) %*% half.D)
#     if(normalised == FALSE) return(D - A)
# 
#   }
