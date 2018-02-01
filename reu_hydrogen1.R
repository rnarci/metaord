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


size_fraction = "0.8-5"
import_data(size_fraction, samples = NULL)

############################################################ Clustering

############################## MDS + k-means

###### Abundance matrices

l = 5

ARI1 = matrix(NA,100,50)
RI1 = matrix(NA,100,50)

for(k in 1:50)
{
  fit1 = cmdscale(jaccard_abundance, eig=TRUE, k=k)
  fit2 = cmdscale(ab_sorensen_abundance,eig=TRUE, k=k)

  for(nb_run in 1:100)
  {
    res1 = kmeans(fit1$points,l)$cluster
    res2 = kmeans(fit2$points,l)$cluster
      
    ARI1[nb_run,k] = adjustedRand(get(paste("res",1,sep="")),get(paste("res",2,sep="")),randMethod="Rand")
    RI1[nb_run,k] = rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
  }
}

ARI = c()
RI = c()

for(k in 1:50){
  ARI[k] = max(ARI1[,k]) 
  RI[k] = max(RI1[,k]) 
}




###### Prevalence matrices


ARI1 = matrix(NA,100,50)
RI1 = matrix(NA,100,50)

for(k in 1:50)
{
  fit1 = cmdscale(kulczynski_prevalence,eig=TRUE, k=k)
  fit2 = cmdscale(ochiai_prevalence,eig=TRUE, k=k)
  
  for(nb_run in 1:100)
  {
    res1 = kmeans(fit1$points,l)$cluster
    res2 = kmeans(fit2$points,l)$cluster
    
    ARI1[nb_run,k] = adjustedRand(get(paste("res",1,sep="")),get(paste("res",2,sep="")),randMethod="Rand")
    RI1[nb_run,k] = rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
  }
}

ARI = c()
RI = c()

for(k in 1:50){
  ARI[k] = max(ARI1[,k]) 
  RI[k] = max(RI1[,k]) 
}




############################## k-NNG + spectral clustering


###### Abundance matrices

l = 5

ARI1 = matrix(NA,100,50)
RI1 = matrix(NA,100,50)

for(k in 1:50)
{
  kNN1 =  make.kNNG(jaccard_abundance, k = k, symm = TRUE, weight = FALSE)
  kNN2 = make.kNNG(ab_sorensen_abundance, k = k, symm = TRUE, weight = FALSE)
  
  for(nb_run in 1:100)
  {
    res1 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
    res2 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
    
    ARI1[nb_run,k] = adjustedRand(get(paste("res",1,sep="")),get(paste("res",2,sep="")),randMethod="Rand")
    RI1[nb_run,k] = rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
  }
}

ARI = c()
RI = c()

for(k in 1:50){
  ARI[k] = max(ARI1[,k]) 
  RI[k] = max(RI1[,k]) 
}


###### Prevalence matrices


ARI1 = matrix(NA,100,50)
RI1 = matrix(NA,100,50)

for(k in 1:50)
{
  kNN1 =  make.kNNG(ochiai_prevalence, k = k, symm = TRUE, weight = FALSE)
  kNN2 = make.kNNG(kulczynski_prevalence, k = k, symm = TRUE, weight = FALSE)
  
  for(nb_run in 1:100)
  {
    res1 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
    res2 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
    
    ARI1[nb_run,k] = adjustedRand(get(paste("res",1,sep="")),get(paste("res",2,sep="")),randMethod="Rand")
    RI1[nb_run,k] = rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
  }
}

ARI = c()
RI = c()

for(k in 1:50){
  ARI[k] = max(ARI1[,k]) 
  RI[k] = max(RI1[,k]) 
}




##################################### Comparison between clusterings 

######### Abundance matrices

fit1 = cmdscale(jaccard_abundance, eig=TRUE, k=10)
fit2 = cmdscale(ab_sorensen_abundance,eig=TRUE, k=10)
res1 = kmeans(fit1$points,l)$cluster
res2 = kmeans(fit2$points,l)$cluster

kNN1 =  make.kNNG(jaccard_abundance, k = 41, symm = TRUE, weight = FALSE)
kNN2 = make.kNNG(ab_sorensen_abundance, k = 41, symm = TRUE, weight = FALSE)
res3 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res4 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
rand.index(get(paste("res",3,sep="")),get(paste("res",4,sep="")))
rand.index(get(paste("res",1,sep="")),get(paste("res",3,sep="")))
rand.index(get(paste("res",2,sep="")),get(paste("res",4,sep="")))

######### Prevalence matrices

fit1 = cmdscale(ochiai_prevalence, eig=TRUE, k=10)
fit2 = cmdscale(kulczynski_prevalence,eig=TRUE, k=10)
res1 = kmeans(fit1$points,l)$cluster
res2 = kmeans(fit2$points,l)$cluster

kNN1 =  make.kNNG(ochiai_prevalence, k = 35, symm = TRUE, weight = FALSE)
kNN2 = make.kNNG(kulczynski_prevalence, k = 35, symm = TRUE, weight = FALSE)
res3 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res4 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
rand.index(get(paste("res",3,sep="")),get(paste("res",4,sep="")))
rand.index(get(paste("res",1,sep="")),get(paste("res",3,sep="")))
rand.index(get(paste("res",2,sep="")),get(paste("res",4,sep="")))

# fit1 = cmdscale(ochiai_prevalence, eig=TRUE, k=10)
# fit2 = cmdscale(ochiai_prevalence,eig=TRUE, k=4)
# fit3 = cmdscale(kulczynski_prevalence, eig=TRUE, k=10)
# fit4 = cmdscale(kulczynski_prevalence,eig=TRUE, k=4)
# res1 = kmeans(fit1$points,l)$cluster
# res2 = kmeans(fit2$points,l)$cluster
# res3 = kmeans(fit3$points,l)$cluster
# res4 = kmeans(fit4$points,l)$cluster
# rand.index(get(paste("res",1,sep="")),get(paste("res",3,sep="")))
# rand.index(get(paste("res",2,sep="")),get(paste("res",4,sep="")))
# rand.index(get(paste("res",1,sep="")),get(paste("res",2,sep="")))
# rand.index(get(paste("res",3,sep="")),get(paste("res",4,sep="")))


############################### To search which k maximises similarity between graphs

############# Abundance matrices

res = c()

for(k in 1:50){
  
  kNN1 = make.kNNG(jaccard_abundance, k = k, symm = TRUE, weight = FALSE)

  kNN2 = make.kNNG(ab_jaccard_abundance, k = k, symm = TRUE, weight = FALSE)

  kNN3 = make.kNNG(braycurtis_abundance, k = k, symm = TRUE, weight = FALSE)

  kNN4 = make.kNNG(ab_ochiai_abundance, k = k, symm = TRUE, weight = FALSE)

  kNN5 = make.kNNG(ab_sorensen_abundance, k = k, symm = TRUE, weight = FALSE)

  kNN6 = make.kNNG(simka_jaccard_abundance, k = k, symm = TRUE, weight = FALSE)
  
  index = 0
  
  for(j in 1:5){
    for(m in (j+1):6){
      index = index + GARI(get(paste("kNN",j,sep="")),get(paste("kNN",m,sep="")))
    }
  }
  res[k] = index/15
}

res
which.max(res)
max(res)

############# Prevalence matrices

res = c()

for(k in 1:50){
  
  kNN7 = make.kNNG(chord_prevalence, k = k, symm = TRUE, weight = FALSE)
  
  kNN8 = make.kNNG(jaccard_prevalence, k = k, symm = TRUE, weight = FALSE)
  
  kNN9 = make.kNNG(kulczynski_prevalence, k = k, symm = TRUE, weight = FALSE)
  
  kNN10 = make.kNNG(ochiai_prevalence, k = k, symm = TRUE, weight = FALSE)
  
  kNN11 = make.kNNG(whittaker_prevalence, k = k, symm = TRUE, weight = FALSE)
  
  kNN12 = make.kNNG(simka_jaccard_prevalence, k = k, symm = TRUE, weight = FALSE)
  
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

####################################################################### Multivariate analysis 

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

library(mice)
mice = complete(mice(design,method="norm.predict",m=1))
design$Phosphates = mice$Phosphates
design$NO2NO3 = mice$NO2NO3

######################################## To calculate predicted distances

D = jaccard_abundance
n = dim(D)[1]

parms1 = c(design$Temperature,design$SSD,design$SI_temperature,design$Depth,design$Silicate,design$SI_nitrates)
X1 = matrix(parms1,nrow=n,ncol=6)
parms2 = design$Temperature
X2 = matrix(parms2,nrow=n,ncol=1)

H1 = X1%*%solve(t(X1)%*%X1)%*%t(X1)
H2 = X2%*%solve(t(X2)%*%X2)%*%t(X2)

D_without_MDS1 = matrix(NA,n,n)
D_without_MDS2 = matrix(NA,n,n)

J = diag(rep(1,n)) - matrix(1,n,n)/n
G = -0.5*J%*%D^2%*%J

for(i in 1:n){
  for(j in 1:n){
    D_without_MDS1[i,j] = sqrt(H1[i,]%*%G%*%H1[,i] + H1[j,]%*%G%*%H1[,j] - 2*H1[i,]%*%G%*%H1[,j])
    D_without_MDS2[i,j] = sqrt(H2[i,]%*%G%*%H2[,i] + H2[j,]%*%G%*%H2[,j] - 2*H2[i,]%*%G%*%H2[,j])
  }
}

fit1 = cmdscale(D, eig=TRUE, k=10)
fit2 = cmdscale(D_without_MDS1, eig=TRUE, k=10)
fit3 = cmdscale(D_without_MDS2, eig=TRUE, k=10)

library(loe)
library(fcd)
library(clues)

kNN1 = make.kNNG(D, k = 37, symm = TRUE, weight = FALSE)
kNN2 = make.kNNG(D_without_MDS1, k = 37, symm = TRUE, weight = FALSE)
kNN3 = make.kNNG(D_without_MDS2, k = 37, symm = TRUE, weight = FALSE)

l = 3

res1 = kmeans(fit1$points,l)$cluster
res2 = kmeans(fit2$points,l)$cluster
res3 = kmeans(fit3$points,l)$cluster

res4 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res5 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res6 = spectral.clustering(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)


########## Boxplots : clustering on numeric distances

temperature_data = data.frame(clus.raw = res1,
                              clus.predicted.all.variables = res2,
                              clus.predicted.only.temperature = res3) %>% 
  mutate(Temperature = design$Temperature)

tdata <- gather(temperature_data, key = "Clustering", value = "Group", clus.raw:clus.predicted.only.temperature)

ggplot(tdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 3)

########## Boxplots : clustering on ordinal distances

temperature_data = data.frame(clus.raw = res4,
                              clus.predicted.all.variables = res5,
                              clus.predicted.only.temperature = res6) %>% 
  mutate(Temperature = design$Temperature)

tdata <- gather(temperature_data, key = "Clustering", value = "Group", clus.raw:clus.predicted.only.temperature)

ggplot(tdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 3)

################################## Map resulting clusters on raw distances (ordinal setting)

library(rworldmap)
newmap <- getMap(resolution = "li")

color = c()
labels = c()
cl = res4

for(i in 1:length(cl)){
  if(cl[i]==1){
    color[i] = "blue"
  }
  if(cl[i]==2){
    color[i] = "gold4"
  }
  if(cl[i]==3){
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
text(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], labels = labels[DCM_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)
text(gps$Mean_longitude[SUR_only_indices], gps$Mean_latitude[SUR_only_indices] + 2, labels = labels[SUR_only_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)



########## Map of resulting clusters : raw vs predict. with Temperature only

color = c()
labels = c()
cl = res6

for(i in 1:length(cl)){
  if(cl[i]==1){
    color[i] = "blue"
  }
  if(cl[i]==2){
    color[i] = "red"
  }
  if(cl[i]==3){
    color[i] = "gold4"
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
text(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], labels = labels[DCM_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)
text(gps$Mean_longitude[SUR_only_indices], gps$Mean_latitude[SUR_only_indices] + 2, labels = labels[SUR_only_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)




########## Map of resulting clusters : raw vs predict. with all covariates

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
  if(cl[i]==3){
    color[i] = "gold4"
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
text(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], labels = labels[DCM_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)
text(gps$Mean_longitude[SUR_only_indices], gps$Mean_latitude[SUR_only_indices] + 2, labels = labels[SUR_only_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)

########## Plot a distance matrix as a heatmap

library(ggplot2)
library(reshape2)

plot_dist_as_heatmap <- function(dist, order = NULL, title = NULL,
                                 low = "#B1F756", high = "#132B13") {
  data <- melt(as(dist, "matrix"))
  colnames(data) <- c("x", "y", "distance")
  if (!is.null(order)) {
    data$x <- factor(data$x, levels = order)
    data$y <- factor(data$y, levels = order)
  }
  p <- ggplot(data, aes(x = x, y = y, fill = distance)) + geom_tile()
  p <- p + theme(axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.text.y = element_blank()) +
    scale_fill_gradient(low = low, high = high)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

plot_dist_as_heatmap(as.dist(D), order = order(as.numeric(labels)), title = "Heatmap for raw distance matrix")
plot_dist_as_heatmap(as.dist(D_without_MDS2), order = order(as.numeric(labels)), title = "Heatmap for predicted distance matrix (only temperature)")
plot_dist_as_heatmap(as.dist(D_without_MDS1), order = order(as.numeric(labels)), title = "Heatmap for predicted distance matrix (all variables)")

######################################## To calculate corrected distances

D_without_MDS1 = matrix(NA,n,n)
D_without_MDS2 = matrix(NA,n,n)
I = diag(rep(1,n))

for(i in 1:n){
  for(j in 1:n){
    D_without_MDS1[i,j] = sqrt((I-H1)[i,]%*%G%*%(I-H1)[,i] + (I-H1)[j,]%*%G%*%(I-H1)[,j] - 2*(I-H1)[i,]%*%G%*%(I-H1)[,j])
    D_without_MDS2[i,j] = sqrt((I-H2)[i,]%*%G%*%(I-H2)[,i] + (I-H2)[j,]%*%G%*%(I-H2)[,j] - 2*(I-H2)[i,]%*%G%*%(I-H2)[,j])
  }
}

fit1 = cmdscale(D, eig=TRUE, k=10)
fit2 = cmdscale(D_without_MDS1, eig=TRUE, k=10)
fit3 = cmdscale(D_without_MDS2, eig=TRUE, k=10)

kNN1 = make.kNNG(D, k = 37, symm = TRUE, weight = FALSE)
kNN2 = make.kNNG(D_without_MDS1, k = 37, symm = TRUE, weight = FALSE)
kNN3 = make.kNNG(D_without_MDS2, k = 37, symm = TRUE, weight = FALSE)

res1 = kmeans(fit1$points,l)$cluster
res2 = kmeans(fit2$points,l)$cluster
res3 = kmeans(fit3$points,l)$cluster


res4 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res5 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res6 = spectral.clustering(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)


########## Boxplots : clustering on numeric distances

temperature_data = data.frame(clus.raw = res1,
                              clus.corrected.all.variables = res2,
                              clus.corrected.only.temperature = res3) %>% 
  mutate(Temperature = design$Temperature)

tdata <- gather(temperature_data, key = "Clustering", value = "Group", clus.raw:clus.corrected.only.temperature)

ggplot(tdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 3)


########## Boxplots : clustering on ordinal distances

temperature_data = data.frame(clus.raw = res4,
                              clus.corrected.all.variables = res5,
                              clus.corrected.only.temperature = res6) %>% 
  mutate(Temperature = design$Temperature)

tdata <- gather(temperature_data, key = "Clustering", value = "Group", clus.raw:clus.corrected.only.temperature)

ggplot(tdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 3)


########## Map of resulting clusters : raw vs correct. with Temperature only

color = c()
labels = c()
cl = res6

for(i in 1:length(cl)){
  if(cl[i]==1){
    color[i] = "gold4"
  }
  if(cl[i]==2){
    color[i] = "red"
  }
  if(cl[i]==3){
    color[i] = "blue"
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
text(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], labels = labels[DCM_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)
text(gps$Mean_longitude[SUR_only_indices], gps$Mean_latitude[SUR_only_indices] + 2, labels = labels[SUR_only_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)




########## Map of resulting clusters : raw vs correct. with all covariates

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
  if(cl[i]==3){
    color[i] = "gold4"
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
text(gps$Mean_longitude[DCM_indices], gps$Mean_latitude[DCM_indices], labels = labels[DCM_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)
text(gps$Mean_longitude[SUR_only_indices], gps$Mean_latitude[SUR_only_indices] + 2, labels = labels[SUR_only_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)

########## Plot a distance matrix as a heatmap

plot_dist_as_heatmap(as.dist(D_without_MDS2), order = order(as.numeric(labels)), title = "Heatmap for corrected distance matrix (only temperature)")
plot_dist_as_heatmap(as.dist(D_without_MDS1), order = order(as.numeric(labels)), title = "Heatmap for corrected distance matrix (all variables)")
