rm(list=objects())

library(vegan)

############################################################ Source custom scripts

source('~/metaord/utils.R')

############################################################ Design

data.wd <- "/home/rnarci/Bureau/CDDrnarci/Donnees/"
design = read.table(file=file.path(data.wd, "param_bioadvection.csv"),sep="",header=TRUE)
rownames(design) <- design$Sample

############################################################ Import data

size_fraction = "20-180"
import_data(size_fraction, samples = rownames(design))

############################################################ Subset design

design <- design[metagenomic_sample, ]
design <- design[,-c(1,5,6,7,13)]

############################################################ ADONIS

adonis(as.dist(jaccard_abundance) ~ Depth, data = design)
adonis(as.dist(jaccard_abundance) ~ Temperature, data = design)
adonis(as.dist(jaccard_abundance) ~ Chlorophyll, data = design)
# adonis(as.dist(jaccard_abundance) ~ Phosphates, data = design)
# adonis(as.dist(jaccard_abundance) ~ Iron, data = design)
# adonis(as.dist(jaccard_abundance) ~ Light, data = design)
adonis(as.dist(jaccard_abundance) ~ Silicate, data = design)
adonis(as.dist(jaccard_abundance) ~ SI_temperature, data = design)
adonis(as.dist(jaccard_abundance) ~ SI_nitrates, data = design)
adonis(as.dist(jaccard_abundance) ~ NH4, data = design)
adonis(as.dist(jaccard_abundance) ~ NP, data = design)
# adonis(as.dist(jaccard_abundance) ~ NO2NO3, data = design)
adonis(as.dist(jaccard_abundance) ~ SSD, data = design)


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
