## Joe Gage
## 21 September 2016

## Run GWAS on slope and MSE values for G x E data

set.seed(1869)
source("http://www.bioconductor.org/biocLite.R")
biocLite("multtest", suppressUpdates=TRUE)

library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library(scatterplot3d)

##source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("./gapit_offline.R")
##source("http://zzlab.net/GAPIT/emma.txt")
source("emma_offline.R")

mse <- read.csv("mse_hybrids_forGWAS.csv", stringsAsFactors=FALSE)
slope <- read.csv("slope_hybrids_forGWAS.csv", stringsAsFactors=FALSE)

if(all(mse[,1] == slope[,1])){
    colnames(mse)[2:ncol(mse)] <- paste0(colnames(mse)[2:ncol(mse)], "_mse")
    colnames(slope)[2:ncol(slope)] <- paste0(colnames(slope)[2:ncol(slope)], "_slope")
    Y <- cbind(slope, mse[,2:ncol(mse)])
}

load("hybrid_geno_forGAPIT.rda")
#K <- read.csv("Kin.VanRaden.csv", stringsAsFactors=FALSE, header=FALSE)

subset=TRUE

## Compression groups chosen in initial model-selection run
##compression <- c(63, 62, 76, 91, 129, 103, 102, 46, 41, 179)
compression <- c(53, 62, 96, 91, 129, 93, 102, 46, 311, 49)
names(compression) <- colnames(Y)[2:ncol(Y)]


for(trait in 2:ncol(Y)){
    if(subset == TRUE){
        keep = Y[!is.na(Y[,trait]),1]
    } else {
        keep = Y[,1]
    }

    groups <- compression[trait-1]

    dir <- paste0("~/Dropbox/papers/GxE/analysis/snp_class/vanRaden/gapit_", colnames(Y)[trait], "_defGroups")
        dir.create(dir)
        setwd(dir)
        gwas <- GAPIT(Y=Y[Y[,1]%in%keep,c(1,trait)],
                      GD=GD[GD[,1]%in%keep, ],
                      GM=GM,
                      PCA.total=0,
                      SNP.fraction=1,
		      SNP.MAF=.005,
                      group.from=groups,
                      group.to=groups
                      )
}
