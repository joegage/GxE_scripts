## Joe Gage
## 17 June 2015

## Consolidate regression parameters from pheno_prep.R into fewer files
## to facilitate GWAS

source("~/Dropbox/Lab/R_functions/alnum.R")

folders <- list.dirs(full.names=T)
files <- list.files(path=folders, pattern="regression_.*\\.txt", full.names=T)

## Read and format raw data
raw <- read.csv("./files/hybrid_data_clean_060915.csv", header=T, stringsAsFactors=F)
traitCols <- c("Experiment", "Rep", "Pedigree", "Pollen.DAP..days.", "Silk.DAP..days.", "Plant.height..in.", "Ear.height..in.", "Grain.yield..bu.A.")
raw <- raw[,traitCols]
colnames(raw) <- c("experiment", "rep", "pedigree", "pollendap", "silkdap", "pltht", "earht", "yld_bu_ac")

index = 1
for(file in files){
    dat <- read.table(file, header=T, stringsAsFactors=F)
    print(dim(dat))

    if(file == files[1]){
        intercept <- matrix(nrow=nrow(dat), ncol=length(files))
        slope <- matrix(nrow=nrow(dat), ncol=length(files))
        mse <- matrix(nrow=nrow(dat), ncol=length(files))
        Rsq <- matrix(nrow=nrow(dat), ncol=length(files))
    }

    ## Drop regression parameters for hybrids planted at less than
    # 3 locations or with less than 6 data points

    remove <- which(dat$nenv < 4 | dat$nobs < 6)
    dat[remove, c("beta0", "beta1", "MSE", "Rsq", "nenv", "nobs")] <- NA

    trait <- gsub(".*/regression_(.*)\\.txt","\\1", file)

    ## Plot distributions of cleaned data
    pdf(paste0("./", trait, "/stability_histograms_", trait, "_clean.pdf"), width=6, height=6)
    par(mfrow=c(3,2))
    hist(dat[,3], breaks=20, col="aquamarine4", main="Slope", ylab="Count", xlab="")
    hist(dat[,4], breaks=20, col="aquamarine4", main="MSE", ylab="Count", xlab="")
    hist(dat[,7], breaks=20, col="aquamarine4", main="R squared", ylab="Count", xlab="")
    hist(dat[,5], breaks=20, col="aquamarine4", main="Number of Obs", ylab="Count", xlab="")
    hist(dat[,6], breaks=20, col="aquamarine4", main="Number of Environments", ylab="Count", xlab="")
    plot(dat[,3], dat[,4], pch=19, ylab="MSE", xlab="Slope")
    dev.off()

    ## Plot regression lines of cleaned data
    pdf(paste0("./", trait, "/stability_regressions_", trait, "_clean.pdf"), width=8, height=4)
    spread <- range(raw[,trait], na.rm=T)
    plot(1, xlim=c((min(spread)*.9),(max(spread)*1.1)), ylim=c((min(spread)*.9),(max(spread)*1.1)),
         xlab="Environmental Mean",
         ylab="Hybrid Phenotype",
         cex.lab=1.5)
    for(i in 1:nrow(dat)){
        if(!is.na(dat[i,"beta1"])){
            abline(a=dat[i,"beta0"], b=(dat[i,"beta1"]+1)) ## +1 for check correction
        }}
    abline(a=0, b=1, col="red")
    dev.off()

    slope[, index] <- dat[,"beta1"]
    intercept[, index] <- dat[,"beta0"]
    mse[, index] <- dat[,"MSE"]
    Rsq[,index] <- dat[,"Rsq"]

    index = index + 1
}

traits <- gsub(".*/regression_(.*)\\.txt","\\1", files)
colnames(slope) <- traits; rownames(slope) = tolower(alnum(dat[,1]))
colnames(intercept) <- traits; rownames(intercept) = tolower(alnum(dat[,1]))
colnames(mse) <- traits; rownames(mse) = tolower(alnum(dat[,1]))
colnames(Rsq) <- traits; rownames(Rsq) = tolower(alnum(dat[,1]))

write.table(slope, "slope_hybrids_forGWAS.csv", row.names=T, col.names=NA, sep=",")
write.table(mse, "mse_hybrids_forGWAS.csv", row.names=T, col.names=NA, sep=",")
