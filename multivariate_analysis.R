rm(list=objects())

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)

############################################################ Source custom scripts

source('~/metaord/utils.R')

############################################################ Design

data.wd <- "/home/rnarci/Bureau/CDDrnarci/Donnees/"
design = read.table(file=file.path(data.wd, "param_bioadvection.csv"),sep="",header=TRUE)
rownames(design) <- design$Sample

############################################################ Import data

size_fraction = "0.8-5"
import_data(size_fraction, samples = rownames(design))

############################################################ Subset design

design <- design[metagenomic_sample, ]
design <- design[,-c(1,6,7)]
# design <- design[,-c(1,5,6,7,13)]

library(mice)
mice = complete(mice(design,method="norm.predict",m=1))
design$Phosphates = mice$Phosphates
design$NO2NO3 = mice$NO2NO3

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
  
  return(c("frac.neg" = mean(eigvalues < 0), ## fraction of negative eigenvalues 
           summary(eigvalues) ## summary of eigenvalues
           ))
}

sapply(matrices.list2, eigen_summary) %>% t()
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

####### Correction euclidienne

D = jaccard_abundance
n = dim(D)[1]

### Par MDS : 

D_MDS1 = matrix(NA,n,n)
D_MDS2 = matrix(NA,n,n)

mod1 = capscale(as.dist(D) ~ Temperature + SSD + SI_temperature + Depth + Silicate + SI_nitrates, data = design)
mod2 = capscale(as.dist(D) ~ Temperature, data = design)
Y_tilde1 = as.matrix(mod1$CCA$u)
Y_tilde1 = scale(Y_tilde1)
Y_tilde2 = as.matrix(mod2$CCA$u)
Y_tilde2 = scale(Y_tilde2)

parms1 = c(design$Temperature,design$SSD,design$SI_temperature,design$Depth,design$Silicate,design$SI_nitrates)
X1 = matrix(parms1,nrow=n,ncol=6)
# parms2 = c(design$Temperature,design$Depth,design$Chlorophyll,design$Phosphates,design$Silicate,design$SI_temperature,design$SI_nitrates,design$NH4,design$NP,design$NO2NO3,design$SSD)
parms2 = design$Temperature
X2 = matrix(parms2,nrow=n,ncol=1)
H1 = X1%*%solve(t(X1)%*%X1)%*%t(X1)
H2 = X2%*%solve(t(X2)%*%X2)%*%t(X2)

Y_hat1 = H1%*%Y_tilde1 
Y_hat2 = H2%*%Y_tilde2 

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

fit1 = cmdscale(D, eig=TRUE, k=4)
fit2 = cmdscale(D_MDS1, eig=TRUE, k=4)
fit3 = cmdscale(D_without_MDS1, eig=TRUE, k=4)
fit4 = cmdscale(D_MDS2, eig=TRUE, k=4)
fit5 = cmdscale(D_without_MDS2, eig=TRUE, k=4)

####### Pre-requis pour clustering ordinal

library(loe)
library(fcd)
library(clues)

kNN1 = make.kNNG(D, k = 35, symm = TRUE, weight = FALSE)
kNN2 = make.kNNG(D_MDS1, k = 35, symm = TRUE, weight = FALSE)
kNN3 = make.kNNG(D_without_MDS1, k = 35, symm = TRUE, weight = FALSE)
kNN4 = make.kNNG(D_MDS2, k = 35, symm = TRUE, weight = FALSE)
kNN5 = make.kNNG(D_without_MDS2, k = 35, symm = TRUE, weight = FALSE)

####### Clustering de l clusters

l = 3

res1 = kmeans(fit1$points,l)$cluster
res2 = kmeans(fit2$points,l)$cluster
res3 = kmeans(fit3$points,l)$cluster
res4 = kmeans(fit4$points,l)$cluster
res5 = kmeans(fit5$points,l)$cluster

res6 = spectral.clustering(kNN1, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res7 = spectral.clustering(kNN2, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res8 = spectral.clustering(kNN3, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res9 = spectral.clustering(kNN4, normalised = TRUE, score = FALSE, K = l, adj = FALSE)
res10 = spectral.clustering(kNN5, normalised = TRUE, score = FALSE, K = l, adj = FALSE)

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

################################## Classical clustering

###### Raw data vs corrected data 1 (only Temperature)

temperature_data = data.frame(clus.raw = res1,
                              clus.corrected.only.temperature = res2,
                              clus.corrected.all.variables = res3) %>% 
  mutate(Temperature = design$Temperature)
head(temperature_data)

tdata <- gather(temperature_data, key = "Clustering", value = "Group", clus.raw:clus.corrected.all.variables)
head(tdata)

ggplot(tdata, aes(x = Group, y = Temperature, group = Group)) + geom_boxplot() + facet_wrap(~Clustering, ncol = 3)

