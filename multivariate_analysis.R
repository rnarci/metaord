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

# missDummy <- function(t)
# {
#   x <- dim(length(t))
#   x[which(!is.na(t))] = 1
#   x[which(is.na(t))] = 0
#   return(x)
# }
# y = as.vector(design$Phosphates)
# ry = missDummy(design$Phosphates)
# x = as.matrix(design)
# mice = mice.impute.norm.predict(y=y,ry=ry,x=x,wy=1-ry)

############################################################ ADONIS

adonis(as.dist(jaccard_abundance) ~ Latitude, data = design)
adonis(as.dist(jaccard_abundance) ~ Longitude, data = design)
adonis(as.dist(jaccard_abundance) ~ Depth, data = design)
adonis(as.dist(jaccard_abundance) ~ Temperature, data = design)
adonis(as.dist(jaccard_abundance) ~ Chlorophyll, data = design)
adonis(as.dist(jaccard_abundance) ~ Phosphates, data = design)
# adonis(as.dist(jaccard_abundance) ~ Iron, data = design)
# adonis(as.dist(jaccard_abundance) ~ Light, data = design)
adonis(as.dist(jaccard_abundance) ~ Silicate, data = design)
adonis(as.dist(jaccard_abundance) ~ SI_temperature, data = design) # SI = seasonality index
adonis(as.dist(jaccard_abundance) ~ SI_nitrates, data = design)
adonis(as.dist(jaccard_abundance) ~ NH4, data = design) # ammonium
adonis(as.dist(jaccard_abundance) ~ NP, data = design) # nitrate/phosphate
adonis(as.dist(jaccard_abundance) ~ NO2NO3, data = design) # nitrite/nitrates
adonis(as.dist(jaccard_abundance) ~ SSD, data = design) # salinity ?


adonis(as.dist(jaccard_abundance) ~ Depth + Temperature + Chlorophyll + Silicate + SI_temperature + SI_nitrates + NH4 + NP + SSD, data = design)
anova.cca(capscale(as.dist(jaccard_abundance) ~ Depth + Temperature + Chlorophyll + Silicate + SI_temperature + SI_nitrates + NH4 + NP + SSD, data = design),by="terms")

##### 0.8-5

adonis(as.dist(whittaker_prevalence) ~ Temperature + SSD + SI_temperature + Depth + Silicate + SI_nitrates + Temperature:SI_nitrates + Temperature:SSD, data = design)
adonis(as.dist(jaccard_abundance) ~ Temperature + SSD + SI_temperature + Depth + Silicate + SI_nitrates + Temperature:SI_nitrates + Temperature:SSD, data = design)
anova.cca(capscale(as.dist(jaccard_abundance) ~ Temperature + SSD + SI_temperature + Depth + Silicate + SI_nitrates + Temperature:SI_nitrates + Temperature:SSD, data = design),by="terms")

##### 5-20

adonis(as.dist(whittaker_prevalence) ~ Temperature + SSD + SI_nitrates + Silicate + SSD:SI_nitrates + Temperature:Silicate, data = design)
adonis(as.dist(jaccard_abundance) ~ Temperature + SSD + SI_nitrates + Silicate + SSD:SI_nitrates + Temperature:Silicate, data = design)
anova.cca(capscale(as.dist(jaccard_abundance) ~ Temperature + SSD + SI_nitrates + Silicate + SSD:SI_nitrates + Temperature:Silicate, data = design),by="terms")

##### 180-2000

adonis(as.dist(whittaker_prevalence) ~ Temperature  + Depth + SSD + Silicate + SI_temperature + NH4 + Temperature:SSD + Temperature:Silicate + NH4:SSD, data = design)
adonis(as.dist(jaccard_abundance) ~ Temperature  + Depth + SSD + Silicate + SI_temperature + NH4 + Temperature:SSD + Temperature:Silicate + NH4:SSD, data = design)
anova.cca(capscale(as.dist(jaccard_abundance) ~ Temperature  + Depth + SSD + Silicate + SI_temperature + NH4 + Temperature:SSD + Temperature:Silicate + NH4:SSD, data = design),by="terms")


##### 0.22-3

adonis(as.dist(whittaker_prevalence) ~ Temperature + SSD + Depth + Chlorophyll + Silicate + NH4 + Temperature:NH4, data = design) # BEST ? 
adonis(as.dist(jaccard_abundance) ~ Temperature + SSD + Depth + Chlorophyll + Silicate + NH4 + Temperature:NH4, data = design)
anova.cca(capscale(as.dist(jaccard_abundance) ~ Temperature + SSD + Depth + Chlorophyll + Silicate + NH4 + Temperature:NH4, data = design),by="terms")

##### 0-0.2

adonis(as.dist(whittaker_prevalence) ~ Temperature + SSD + Depth + NH4 + Silicate + Temperature:SSD + Temperature:Silicate + Temperature:NH4, data = design)
adonis(as.dist(jaccard_abundance) ~ Temperature + SSD + Depth + NH4 + Silicate + Temperature:SSD + Temperature:Silicate + Temperature:NH4, data = design)
anova.cca(capscale(as.dist(jaccard_abundance) ~ Temperature + SSD + Depth + NH4 + Silicate + Temperature:SSD + Temperature:Silicate + Temperature:NH4, data = design),by="terms")

##### 20-180

adonis(as.dist(whittaker_prevalence) ~ Temperature + Silicate + NH4 + Depth + SI_temperature + Temperature:NH4 + Temperature:Silicate, data = design)
adonis(as.dist(jaccard_abundance) ~ Temperature + Silicate + NH4 + Depth + SI_temperature + Temperature:NH4 + Temperature:Silicate, data = design)
anova.cca(capscale(as.dist(jaccard_abundance) ~ Temperature + Silicate + NH4 + Depth + SI_temperature + Temperature:NH4 + Temperature:Silicate, data = design),by="terms")



########### Plot des stations en fonction de la longitude et la latitude



########## Points de la reunion du 09/01/2018 : focalisation sur la fraction de taille 0.8-5

## Eloignement de la geometrie euclidienne ? 

matrices.list1 <- c("jaccard_abundance", "ab_jaccard_abundance", "braycurtis_abundance", 
                   "ab_ochiai_abundance", "ab_sorensen_abundance", "simka_jaccard_abundance", 
                   "chord_prevalence", "jaccard_prevalence", "kulczynski_prevalence", 
                   "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence")

matrices.list2 <- c("jaccard_abundance", "ochiai_abundance", "sorensen_abundance", 
                   "simka_jaccard_abundance", "chord_hellinger_prevalence", 
                   "jaccard_canberra_prevalence", "kulczynski_prevalence", 
                   "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence",
                   "sorensen_braycurtis_prevalence")

# n = dim(jaccard_abundance)[1]
# H = diag(rep(1,n)) - matrix(1,n,n)/n

eigen_summary <- function(mat) {
  D <- get(mat)
  n = dim(D)[1]
  J = diag(rep(1,n)) - matrix(1,n,n)/n
  G = -0.5*J%*%D^2%*%J
  eigvalues <- eigen(G, only.values = T)$values
  
  return(c("nb.zero.val" = sum(eigvalues >= - 10^-8 & eigvalues <= 10^-8),"nb.pos.val" = sum(eigvalues > 10^-8),"nb.neg.val" = sum(eigvalues < - 10^-8), ## number of negative eigenvalues 
           round(summary(eigvalues), digits = 2) ## summary of eigenvalues
           ))
}

sapply(matrices.list1, eigen_summary) %>% t()
## t(sapply(matrices.list1, eigen_summary))

# for(mat in matrices.list1){
# D = get(mat)
# 
# # Id = diag(identity(n))
# # un = t(t(rep(1,n)))
# # scale = Id-(1/n)*un%*%t(un)
# 
# A = -0.5*H%*%D^2%*%H
# print(summary(eigen(A)$values))# = total sum of squares
# }
# 
# D <- matrix(c(0, 2, 2, 1, 
#               2, 0, 2, 1, 
#               2, 2, 0, 1,
#               1, 1, 1, 0), 
#             nrow = 4, ncol = 4)
# n <- nrow(D); H = diag(rep(1,n)) - matrix(1,n,n)/n
# 
# summary(eigen(-1/2 * H %*% D^2 %*% H, only.values = TRUE))
# 
# # jaccard_abunœdance : fraction 0.8-5
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.0000  0.4369  0.4512  0.4741  0.4744  1.2059 

## Etude des liens entre les variables

library(GGally)

ggpairs(design)

# corr > 0.55 : Temperature:SSD ; Temperature:Silicate ; Silicate:SSD ; SI_nitrates:NH4 ; Phosphates:SSD ; Phosphates:NO2NO3
# Phosphates:Silicate ; NO2NO3:SSD ; NO2NO3:Silicate 

# Multi-colinearite (ie det(X'X) proche de 0) ? 

X = matrix(NA,nrow=dim(design)[1],ncol=dim(design)[2])

for(p in 1:dim(design)[2]){
  X[,p] = design[,p]
}

det(t(X)%*%X) # de l'ordre de 10^31

#### Tolerance

library(plm)

tol_lm1 = 1 - (summary(lm(Temperature ~ ., design)))$r.squared
tol_lm2 = 1 - (summary(lm(Depth ~ ., design)))$r.squared
tol_lm3 = 1 - (summary(lm(Chlorophyll ~ ., design)))$r.squared
tol_lm4 = 1 - (summary(lm(Silicate ~ ., design)))$r.squared
tol_lm5 = 1 - (summary(lm(SI_temperature ~ ., design)))$r.squared
tol_lm6 = 1 - (summary(lm(SI_nitrates ~ ., design)))$r.squared
tol_lm7 = 1 - (summary(lm(NH4 ~ ., design)))$r.squared
tol_lm8 = 1 - (summary(lm(Phosphates ~ ., design)))$r.squared
tol_lm9 = 1 - (summary(lm(NO2NO3 ~ ., design)))$r.squared
tol_lm10 = 1 - (summary(lm(SSD ~ ., design)))$r.squared
tol_lm11 = 1 - (summary(lm(NP ~ ., design)))$r.squared

## Transformation de variables

# Premier facon

boxplot(design$Temperature)
boxplot(design$Depth)
boxplot(design$Chlorophyll)
boxplot(design$Silicate)
boxplot(design$SI_temperature)
boxplot(design$SI_nitrates)
boxplot(design$SSD)
boxplot(design$NH4)
boxplot(design$NO2NO3)
boxplot(design$Phosphates)
boxplot(design$NP)


histogram(design$Temperature,xlab="Temperature") ; histogram(design$Temperature^2) ; histogram(1/design$Temperature)
histogram(sqrt(design$Temperature)) ; histogram(log(design$Temperature)) ; histogram(exp(design$Temperature))

histogram(design$Depth,xlab="Depth") ; histogram(design$Depth^2) ; histogram(1/design$Depth)
histogram(sqrt(design$Depth)) ; histogram(log(design$Depth)) ;histogram(exp(design$Depth))

histogram(design$Chlorophyll,xlab="Chlorophyll") ; histogram(design$Chlorophyll^2) ; histogram(1/(design$Chlorophyll+1),xlab="1/(Chlorophyll+1)")
histogram(sqrt(design$Chlorophyll)) ; histogram(log(design$Chlorophyll)) ;histogram(exp(design$Chlorophyll))

histogram(design$Silicate,xlab="Silicate") ; histogram(design$Silicate^2) ; histogram(1/design$Silicate)
histogram(sqrt(design$Silicate)) ; histogram(log(design$Silicate)) ;histogram(exp(design$Silicate))

histogram(design$SI_temperature,xlab="SI_Temperature") ; histogram(design$SI_temperature^2) ; histogram(1/design$SI_temperature)
histogram(sqrt(design$SI_temperature)) ; histogram(log(design$SI_temperature)) ;histogram(exp(design$SI_temperature))
# shapiro.test(log(design$SI_temperature)) = 7 %

histogram(design$SI_nitrates,xlab="SI_Nitrates") ; histogram(design$SI_nitrates^2) ; histogram(1/design$SI_nitrates)
histogram(sqrt(design$SI_nitrates)) ; histogram(log(design$SI_nitrates)) ;histogram(exp(design$SI_nitrates))

histogram(design$SSD,xlab="SSD") ; histogram(design$SSD^2) ; histogram(1/design$SSD)
histogram(sqrt(design$SSD)) ; histogram(log(design$SSD)) ;histogram(exp(design$SSD))

histogram(design$NH4,xlab="NH4") ; histogram(design$NH4^2) ; histogram(1/design$NH4,xlab="1/NH4")
histogram(sqrt(design$NH4)) ; histogram(log(design$NH4)) ;histogram(exp(design$NH4))

histogram(design$NO2NO3,xlab="NO2NO3") ; histogram(design$NO2NO3^2) ; histogram(1/design$NO2NO3)
histogram(sqrt(design$NO2NO3)) ; histogram(log(design$NO2NO3)) ;histogram(exp(design$NO2NO3))

histogram(design$Phosphates,xlab="Phosphates") ; histogram(design$Phosphates^2) ; histogram(1/design$Phosphates)
histogram(sqrt(design$Phosphates)) ; histogram(log(design$Phosphates)) ;histogram(exp(design$Phosphates))

histogram(design$NP,xlab="NP") ; histogram(design$NP^2) ; histogram(1/design$NP)
histogram(sqrt(design$NP)) ; histogram(log(design$NP),xlab="log(NP)") ;histogram(exp(design$NP))

parm = design$NP
inv.parm = 1/parm

adonis(as.dist(jaccard_abundance) ~ parm)
adonis(as.dist(jaccard_abundance) ~ parm^2)
adonis(as.dist(jaccard_abundance) ~ inv.parm)
adonis(as.dist(jaccard_abundance) ~ sqrt(abs(parm)))
adonis(as.dist(jaccard_abundance) ~ log(abs(parm)))
adonis(as.dist(jaccard_abundance) ~ exp(parm))

# design$Chlorophyll : normal -> F.Model = 1.6337 et R2 = 0.0103 ; 1/(x+1) -> F.model =  1.9635 et R2 = 0.0124
# design$NH4 : normal -> F.Model = 1.7679 et R2 = 0.0111 ; 1/x -> F.model =  1.9368 et R2 = 0.0122
# design$NP : normal -> F.Model = 1.0195 et R2 = 0.0065 ; log(x) -> F.model =  1.6458 et R2 = 0.0104
# Sinon changements minimes pour transformations exp, log, sqrt, ^2 et 1/

# Deuxieme facon (avec utilisation de boxcox)

new_design = design

library(EnvStats)

for(i in 1:dim(design)[2]){
  beta = min(design[,i])
  if(beta <= 0){
    new_design[,i] = design[,i] + 1 - beta
  }
  lambda = boxcox(new_design[,i], lambda = c(-5,5), optimize = TRUE, objective.name = "PPCC")$lambda
  new_design[,i] = boxcoxTransform(new_design[,i], lambda = lambda)
}

ggpairs(new_design)

X = matrix(NA,nrow=dim(design)[1],ncol=dim(design)[2])

for(p in 1:dim(design)[2]){
  X[,p] = new_design[,p]
}

det(t(X)%*%X) # de l'ordre de 10^20

histogram(new_design$Temperature)
histogram(new_design$Depth)
histogram(new_design$Chlorophyll)
histogram(new_design$Silicate)
histogram(new_design$SI_temperature)
histogram(new_design$SI_nitrates)
histogram(new_design$NH4)
histogram(new_design$NP)
histogram(new_design$SSD)
histogram(new_design$NO2NO3)
histogram(new_design$Phosphates)

shapiro.test(new_design$Temperature)
shapiro.test(new_design$Depth)
shapiro.test(new_design$Chlorophyll)
shapiro.test(new_design$Silicate)
shapiro.test(new_design$SI_temperature)
shapiro.test(new_design$SI_nitrates)
shapiro.test(new_design$NH4)
shapiro.test(new_design$NP)
shapiro.test(new_design$SSD)
shapiro.test(new_design$NO2NO3)
shapiro.test(new_design$Phosphates)

adonis(as.dist(jaccard_abundance) ~ Temperature, data = design)
adonis(as.dist(jaccard_abundance) ~ Temperature, data = new_design)
adonis(as.dist(jaccard_abundance) ~ Depth, data = design)
adonis(as.dist(jaccard_abundance) ~ Depth, data = new_design)
adonis(as.dist(jaccard_abundance) ~ Chlorophyll, data = design)
adonis(as.dist(jaccard_abundance) ~ Chlorophyll, data = new_design)
adonis(as.dist(jaccard_abundance) ~ Silicate, data = design)
adonis(as.dist(jaccard_abundance) ~ Silicate, data = new_design)
adonis(as.dist(jaccard_abundance) ~ SI_temperature, data = design)
adonis(as.dist(jaccard_abundance) ~ SI_temperature, data = new_design)
adonis(as.dist(jaccard_abundance) ~ SI_nitrates, data = design)
adonis(as.dist(jaccard_abundance) ~ SI_nitrates, data = new_design)
adonis(as.dist(jaccard_abundance) ~ NH4, data = design)
adonis(as.dist(jaccard_abundance) ~ NH4, data = new_design)
adonis(as.dist(jaccard_abundance) ~ SSD, data = design)
adonis(as.dist(jaccard_abundance) ~ SSD, data = new_design)
adonis(as.dist(jaccard_abundance) ~ NP, data = design)
adonis(as.dist(jaccard_abundance) ~ NP, data = new_design)
adonis(as.dist(jaccard_abundance) ~ NO2NO3, data = design)
adonis(as.dist(jaccard_abundance) ~ NO2NO3, data = new_design)
adonis(as.dist(jaccard_abundance) ~ Phosphates, data = design)
adonis(as.dist(jaccard_abundance) ~ Phosphates, data = new_design)

####################### Distances entre valeurs prédites

D = jaccard_abundance
n = dim(D)[1]

### Par MDS : 

D_MDS1 = matrix(NA,n,n)
D_MDS2 = matrix(NA,n,n)

Y_tilde = cmdscale(D, k = n - 1)

# parms1 = c(design$Temperature,design$SSD,design$SI_temperature,design$Depth,design$Silicate,design$SI_nitrates)
# X1 = matrix(parms1,nrow=n,ncol=6)
parms1 = design$Depth
X1 = matrix(parms1,nrow=n,ncol=1)
parms2 = design$Longitude
X2 = matrix(parms2,nrow=n,ncol=1)

H1 = X1%*%solve(t(X1)%*%X1)%*%t(X1)
H2 = X2%*%solve(t(X2)%*%X2)%*%t(X2)

Y_hat1 = H1%*%Y_tilde
Y_hat2 = H2%*%Y_tilde

for(i in 1:n){
  for(j in 1:n){
D_MDS1[i,j] = sqrt(sum((Y_hat1[i,] - Y_hat1[j,])^2)) 
D_MDS2[i,j] = sqrt(sum((Y_hat2[i,] - Y_hat2[j,])^2)) 
  }
}

## Sans MDS

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

####### Pre-requis pour clustering par k-means

fit1 = cmdscale(D, eig=TRUE, k=10)
fit2 = cmdscale(D_without_MDS1, eig=TRUE, k=10)
fit3 = cmdscale(D_without_MDS2, eig=TRUE, k=10)

####### Pre-requis pour clustering ordinal

library(loe)
library(fcd)
library(clues)

kNN1 = make.kNNG(D, k = 37, symm = TRUE, weight = FALSE)
kNN2 = make.kNNG(D_without_MDS1, k = 37, symm = TRUE, weight = FALSE)
kNN3 = make.kNNG(D_without_MDS2, k = 37, symm = TRUE, weight = FALSE)

####### Clustering de l clusters

l = 2

res1 = kmeans(fit1$points,l, nstart = 1000)$cluster
res2 = kmeans(fit2$points,l, nstart = 1000)$cluster
res3 = kmeans(fit3$points,l, nstart = 1000)$cluster

res4 = spectral.clustering.new(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res5 = spectral.clustering.new(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res6 = spectral.clustering.new(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

# ####### Comparaison des clusterings
# 
# tab = matrix(NA,10,10)
# 
# for(i in 1:10){
#   for(j in 1:10){
#     tab[i,j] = adjustedRand(get(paste("res",i,sep="")),get(paste("res",j,sep="")),randMethod="Rand") 
#   }
# }
# tab

################################## Comparaison des clusterings

###### Example : raw data vs predicted data 1 (only Temperature)

temperature_data = data.frame(clus.raw = res1,
                              clus.predicted.all.variables = res2,
                              clus.predicted.only.temperature = res3) %>% 
  mutate(Temperature = design$Temperature)
head(temperature_data)

tdata <- gather(temperature_data, key = "Clustering", value = "Group", clus.raw:clus.predicted.only.temperature)
head(tdata)

ggplot(tdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 3)




depth_data = data.frame(clus.raw = res1,
                              clus.corrected.all.variables = res2,
                              clus.corrected.only.temperature = res3) %>%
  mutate(Depth = design$Depth)

ddata <- gather(depth_data, key = "Clustering", value = "Group", clus.raw:clus.corrected.only.temperature)

ggplot(ddata, aes(x = Group, y = Depth, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 3)






####################### Distances entre résidus (la fonction capscale n'est pas necessaire : on peut utiliser)
####################### la fonction cdmscale avec une dimension egale a n-1

D = jaccard_abundance
n = dim(D)[1]

### Par MDS : 

D_MDS1 = matrix(NA,n,n)
D_MDS2 = matrix(NA,n,n)

Y_tilde = cmdscale(D,k=n-1)

# parms1 = c(design$Temperature,design$SSD,design$SI_temperature,design$Depth,design$Silicate,design$SI_nitrates)
# X1 = matrix(parms1,nrow=n,ncol=6)
parms1 = c(design$Temperature,design$SSD,design$Silicate)
X1 = matrix(parms1,nrow=n,ncol=3)
parms2 = design$Depth
X2 = matrix(parms2,nrow=n,ncol=1)

H1 = X1%*%solve(t(X1)%*%X1)%*%t(X1)
H2 = X2%*%solve(t(X2)%*%X2)%*%t(X2)

Y_hat1 = H1%*%Y_tilde 
Y_hat2 = H2%*%Y_tilde

R1 = Y_tilde - Y_hat1
R2 = Y_tilde - Y_hat2

for(i in 1:n){
  for(j in 1:n){
    D_MDS1[i,j] = sqrt(sum((R1[i,] - R1[j,])^2)) 
    D_MDS2[i,j] = sqrt(sum((R2[i,] - R2[j,])^2)) 
  }
}

## Sans MDS

D_without_MDS1 = matrix(NA,n,n)
D_without_MDS2 = matrix(NA,n,n)

J = diag(rep(1,n)) - matrix(1,n,n)/n
G = -0.5*J%*%D^2%*%J
I = diag(rep(1,n))
                                                                                                                      
for(i in 1:n){
  for(j in 1:n){
    D_without_MDS1[i,j] = sqrt((I-H1)[i,]%*%G%*%(I-H1)[,i] + (I-H1)[j,]%*%G%*%(I-H1)[,j] - 2*(I-H1)[i,]%*%G%*%(I-H1)[,j])
    D_without_MDS2[i,j] = sqrt((I-H2)[i,]%*%G%*%(I-H2)[,i] + (I-H2)[j,]%*%G%*%(I-H2)[,j] - 2*(I-H2)[i,]%*%G%*%(I-H2)[,j])
  }
}

####### Pre-requis pour clustering par k-means

fit1 = cmdscale(D, eig=TRUE, k=10)
fit2 = cmdscale(D_without_MDS1, eig=TRUE, k=10)
fit3 = cmdscale(D_without_MDS2, eig=TRUE, k=10)

####### Pre-requis pour clustering ordinal

library(loe)
library(fcd)
library(clues)

kNN1 = make.kNNG(D, k = 37, symm = TRUE, weight = FALSE)
kNN2 = make.kNNG(D_without_MDS1, k = 37, symm = TRUE, weight = FALSE)
kNN3 = make.kNNG(D_without_MDS2, k = 37, symm = TRUE, weight = FALSE)

####### Clustering de l clusters

l = 2

res1 = kmeans(fit1$points,l, nstart = 1000)$cluster
res2 = kmeans(fit2$points,l, nstart = 1000)$cluster
res3 = kmeans(fit3$points,l, nstart = 1000)$cluster


res4 = spectral.clustering.new(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res5 = spectral.clustering.new(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res6 = spectral.clustering.new(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)


# ####### Comparaison des clusterings
# 
# tab = matrix(NA,10,10)
# 
# for(i in 1:10){
#   for(j in 1:10){
#     tab[i,j] = adjustedRand(get(paste("res",i,sep="")),get(paste("res",j,sep="")),randMethod="Rand") 
#   }
# }
# tab

################################## Comparaison des clusterings


temperature_data = data.frame(clus.raw = res1,
                              clus.corrected.all.variables = res2,
                              clus.corrected.only.temperature = res3) %>% 
  mutate(Temperature = design$Temperature)
head(temperature_data)

tdata <- gather(temperature_data, key = "Clustering", value = "Group", clus.raw:clus.corrected.only.temperature)
head(tdata)

ggplot(tdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 3)

temperature_data = data.frame(clus.raw = res4,
                              clus.corrected.all.variables = res5,
                              clus.corrected.only.temperature = res6) %>% 
  mutate(Temperature = design$Temperature)
head(temperature_data)

tdata <- gather(temperature_data, key = "Clustering", value = "Group", clus.raw:clus.corrected.only.temperature)
head(tdata)

ggplot(tdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 3)

# depth_data = data.frame(clus.raw = res1,
#                               clus.corrected.all.variables = res2,
#                               clus.corrected.only.temperature = res3) %>% 
#   mutate(Depth = design$Depth)
# 
# ddata <- gather(depth_data, key = "Clustering", value = "Group", clus.raw:clus.corrected.only.temperature)
# 
# ggplot(ddata, aes(x = Group, y = Depth, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 3)

################################### Representation graphique : longitude en fonction de latitude

######### Premiere facon

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

######### Deuxieme facon

library(ggmap)
library(maps)
library(mapdata)


sbbox <- make_bbox(lon = gps$Mean_longitude, lat = gps$Mean_latitude, f = .1)
sq_map <- get_map(location = sbbox, maptype = "watercolor", source = "google")

color = c()
label = c()

for(i in 1:length(res4)){
  if(res4[i]==1){
    color[i] = "black"
  }
  if(res4[i]==2){
    color[i] = "red"
  }
  if(res4[i]==3){
    color[i] = "blue"
  }
  if(str_length(rownames(gps)[i]) == 5){
    label[i] = substr(rownames(gps)[i],1,1)
  }
  if(str_length(rownames(gps)[i]) == 6){
    label[i] = substr(rownames(gps)[i],1,2)
  }
  if(str_length(rownames(gps)[i]) == 7){
    label[î] = substr(rownames(gps)[i],1,3)
  }
}

ggmap(sq_map) + geom_point(data = gps, mapping = aes(x = gps$Mean_longitude, y = gps$Mean_latitude), color = color, pch = 17, cex = 2) # +
  #geom_text(data = gps, aes(label = rownames(gps)), angle = 60, hjust = 0, color = "yellow")


######## Plot of two distance matrices

plot(as.dist(D),as.dist(D_without_MDS1),xlab="Distances brutes",ylab="Distances prédites",pch=16,col="blue", main="Toutes les covariables",cex.main=1.8,cex.lab=1.3)
plot(as.dist(D),as.dist(D_without_MDS2),xlab="Distances brutes",ylab="Distances prédites",pch=16,col="blue", main="Temperature uniquement",cex.main=1.8,cex.lab=1.3)
plot(as.dist(D),as.dist(lagrangien),xlab="Distances brutes",ylab="Distances Lagrangiennes",pch=16,col="blue",cex.lab=1.3)

library(ggplot2)
library(reshape2)

## Plot a distance matrix as a heatmap with samples sorted according to
## order vector
plot_dist_as_heatmap <- function(dist, order = NULL, title = NULL,
                                 low = "#B1F756", high = "#132B13") {
  ## Args:
  ## - dist: distance matrix (dist class)
  ## - order: (optional) ordering of the samples of dist for representation
  ## - title: (optional) graph title
  ## - low, high: (optional) Colours for low and high ends of the gradient
  ##
  ## Returns:
  ## - a ggplot2 object
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

